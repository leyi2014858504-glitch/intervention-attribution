"""M2+M3:MCTS(渐进加宽 + 转置表 + PUCT/策略先验 + value-net 截断)。

选择策略两档(可切换,消融必需):
- **UCT** (`use_puct=False`):`Q + c·sqrt(ln N(s)/N(s,a))`
- **PUCT**(`use_puct=True`):`Q + c_puct·π(a|s)·sqrt(N(s))/(1+N(s,a))`

PUCT 的探索项衰减更快(1/N 而非 1/√N),更偏利用 ⇒ 树更容易往深处走;
`π(a|s)` 由 `PolicyPrior` 给出(uniform / heuristic / learned)。

连续动作靠**固定分支宽度**(每节点展开 `width` 个候选,候选数**不随访问数增长**)处理:
`π` 决定**新候选从哪采**,其权重同时进入 PUCT。值一律由箱模型 rollout / value-net 兜底
—— 近似只在先验里。
**已去除渐进加宽**:原 ⌊k·N^α⌋ 随 N 膨胀 ⇒ 预算全被根吃掉、树只长到 2 层(plan_len=2),
实测使 MCTS 反而不如 1 步贪心(见 exp_temporal_structure 的浅树混淆项)。
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .model import EnvModel
from .policy import PolicyPrior
from .value_net import ValueNet


@dataclass(frozen=True)
class MCTSConfig:
    budget: int = 200                 # 搜索迭代数
    #: **固定分支宽度**(无渐进加宽):每节点最多展开 width 个候选子节点后才继续下潜。
    #: width≥2;=1 会退化成单链(每层只有一个子 ⇒ 无比较)。预算给深度的张力:
    #: 树可及深度 ≈ budget/width,故 width 越大越浅,需配套更大 budget。
    width: int = 16
    c: float = 1.4                    # UCT 探索常数
    seed: int = 0
    rollout_max: Optional[int] = None  # 截断 rollout 步数(None=到 horizon)
    # --- M3-2 PUCT + 策略先验 ---
    use_puct: bool = False
    c_puct: float = 1.4               # PUCT 探索常数
    policy_mode: str = "uniform"      # uniform | heuristic | learned
    policy_sigma: float = 0.30
    policy_min_samples: int = 24
    policy_retrain_every: int = 40
    # --- M3-1 value-net 截断 ---
    use_value_net: bool = False
    value_w: float = 1.0              # w:1=纯网络(零 rollout),0=纯 rollout
    warmup_iters: int = 50            # 前 N 次迭代强制完整 rollout 以收集标签
    retrain_every: int = 50           # 每 N 次迭代重训一次
    label_fraction: float = 0.1       # 就绪后仍以该概率做完整 rollout 刷新标签
    vnet_min_samples: int = 32
    vnet_hidden: tuple = (32, 32)


@dataclass
class Node:
    state: np.ndarray
    step_idx: int
    N: int = 0
    W: float = 0.0                   # 累积回报之和
    children: dict[int, "Node"] = field(default_factory=dict)
    actions: dict[int, np.ndarray] = field(default_factory=dict)   # child_id -> 动作
    priors: dict[int, float] = field(default_factory=dict)         # child_id -> π 权重

    @property
    def q(self) -> float:
        return self.W / self.N if self.N else 0.0

    def is_leaf(self) -> bool:
        return len(self.children) == 0


@dataclass
class SearchResult:
    best_action: np.ndarray
    root: Node
    transposition_hits: int
    model_cache: dict
    vnet_samples: int = 0
    vnet_ready: bool = False
    vnet_fits: int = 0
    policy_samples: int = 0
    policy_ready: bool = False
    policy_fits: int = 0


class MCTS:
    def __init__(self, model: EnvModel, config: MCTSConfig | None = None):
        self.model = model
        self.cfg = config or MCTSConfig()
        self.rng = np.random.default_rng(self.cfg.seed)
        self._trans: dict[tuple, Node] = {}
        self.tp_hits = 0
        self.vnet = (ValueNet(model.backend, seed=self.cfg.seed,
                              horizon=model.horizon,
                              min_samples=self.cfg.vnet_min_samples,
                              hidden=self.cfg.vnet_hidden)
                     if self.cfg.use_value_net else None)
        self.policy = (PolicyPrior(model.backend, mode=self.cfg.policy_mode,
                                   sigma=self.cfg.policy_sigma, seed=self.cfg.seed,
                                   min_samples=self.cfg.policy_min_samples,
                                   horizon=model.horizon)
                       if self.cfg.use_puct else None)

    # ---- 转置表键(状态哈希 + 深度) ----
    def _key(self, state: np.ndarray, step_idx: int) -> tuple:
        return (step_idx, np.asarray(state, float).tobytes())

    def _get_or_make(self, state: np.ndarray, step_idx: int) -> Node:
        key = self._key(state, step_idx)
        node = self._trans.get(key)
        if node is not None:
            self.tp_hits += 1
            return node
        node = Node(state=np.asarray(state, float), step_idx=step_idx)
        self._trans[key] = node
        return node

    def _branch_width(self) -> int:
        """固定分支宽度(已去除渐进加宽):候选数不随访问数 N 增长。"""
        return max(2, int(self.cfg.width))

    # ---- 选择:UCT 或 PUCT ----
    def _best_child(self, node: Node) -> Node:
        logN = math.log(node.N + 1e-12)
        sqrtN = math.sqrt(max(node.N, 1))
        prior_sum = sum(node.priors.values()) or 1.0
        best, best_val = None, -np.inf
        for cid, ch in node.children.items():
            q = ch.W / ch.N if ch.N else 0.0
            if self.policy is not None:
                p = node.priors.get(cid, 1.0) / prior_sum
                u = self.cfg.c_puct * p * sqrtN / (1.0 + ch.N)
            else:
                u = self.cfg.c * math.sqrt(logN / (ch.N + 1e-12))
            val = q + u
            if val > best_val:
                best_val, best = val, ch
        return best if best is not None else node

    def _descend(self, root: Node) -> list[Node]:
        """从根下降,返回访问路径(停在可扩展的节点或 horizon 处)。"""
        path, cur = [root], root
        while (cur.step_idx < self.model.horizon
               and len(cur.children) >= self._branch_width()):
            cur = self._best_child(cur)
            if cur is root:          # 无子可择(异常保护)
                break
            path.append(cur)
        return path

    def _expand(self, node: Node) -> tuple[Node, np.ndarray, float, bool]:
        if self.policy is not None:
            action = self.policy.sample(node.state, node.step_idx, self.rng)
            prior = self.policy.weight(node.state, node.step_idx, action)
        else:
            action = self.rng.random(self.model.backend.action_dim)
            prior = 1.0
        nxt, r, done = self.model.step(node.state, node.step_idx, action)
        child = self._get_or_make(nxt, node.step_idx + 1)
        cid = len(node.children)
        node.children[cid] = child
        node.actions[cid] = action
        node.priors[cid] = prior
        return child, action, r, done

    def _backprop(self, path: list[Node], value: float) -> None:
        for nd in path:
            nd.N += 1
            nd.W += value

    def search(self, root_state: np.ndarray) -> SearchResult:
        root = self._get_or_make(root_state, 0)
        for it in range(self.cfg.budget):
            path = self._descend(root)
            node = path[-1]
            if node.step_idx < self.model.horizon:
                child, action, r, done = self._expand(node)
                path.append(child)
                value = r + (0.0 if done else self._leaf_value(child, it))
                if self.policy is not None:
                    self.policy.add_sample(node.state, node.step_idx, action, value)
            else:
                value = 0.0
            self._backprop(path, value)
            if self.vnet is not None and (it + 1) % self.cfg.retrain_every == 0:
                self.vnet.fit()
            if self.policy is not None and (it + 1) % self.cfg.policy_retrain_every == 0:
                self.policy.fit()
        if root.children:
            best_cid = max(root.children, key=lambda cid: root.children[cid].N)
            best_action = root.actions[best_cid]
        else:
            best_action = np.full(self.model.backend.action_dim, 0.5)
        return SearchResult(best_action=np.asarray(best_action, float), root=root,
                            transposition_hits=self.tp_hits,
                            model_cache=self.model.cache_stats(),
                            vnet_samples=(self.vnet.n_samples if self.vnet else 0),
                            vnet_ready=(self.vnet.ready if self.vnet else False),
                            vnet_fits=(self.vnet.n_fits if self.vnet else 0),
                            policy_samples=(len(self.policy._V) if self.policy else 0),
                            policy_ready=(self.policy.ready if self.policy else False),
                            policy_fits=(self.policy.n_fits if self.policy else 0))

    def _rollout(self, child: Node) -> float:
        return self.model.rollout(child.state, child.step_idx, rng=self.rng,
                                  max_steps=self.cfg.rollout_max)

    def _leaf_value(self, child: Node, it: int) -> float:
        """叶节点价值:按配置混合 value-net 与真实 rollout;需要时记录训练标签。

        w 的取值(warmup 未过 / 网络未就绪 / 抽到标签刷新 ⇒ w=0):
        - w=0   : 完整 rollout(同时把 (s,真实回报) 存入训练集)
        - w=1   : 纯网络,零 rollout(单次迭代成本 ≈ 1 步)
        """
        vnet = self.vnet
        if vnet is None:
            return self._rollout(child)
        if self.cfg.value_w <= 0.0:
            # 纯 rollout 模式:每次必收集标签,且不额外消耗随机数(与禁用路径等价)
            true_v = self._rollout(child)
            vnet.add_sample(child.state, child.step_idx, true_v)
            return true_v
        collect = self.rng.random() < self.cfg.label_fraction
        w = 0.0
        if it >= self.cfg.warmup_iters and vnet.ready and not collect:
            w = float(self.cfg.value_w)
        if w >= 1.0:
            return vnet.predict(child.state, child.step_idx)
        true_v = self._rollout(child)
        vnet.add_sample(child.state, child.step_idx, true_v)
        if w <= 0.0:
            return true_v
        return w * vnet.predict(child.state, child.step_idx) + (1.0 - w) * true_v

    def best_path(self, root: Node) -> list[np.ndarray]:
        """沿访问数最大的子节点贪心下降,提取开环动作序列(评估规划质量用)。"""
        seq, cur = [], root
        while cur.children and cur.step_idx < self.model.horizon:
            cid = max(cur.children, key=lambda c: cur.children[c].N)
            seq.append(cur.actions[cid])
            cur = cur.children[cid]
        return seq
