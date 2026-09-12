"""atmos:大气化学干预规划 + 归因研究框架(方法验证型)。

后端可替换:MockBackend(测试替身)→ CB05 → SAPRC-07,均实现同一
`ChemistryBackend` 契约。上层(MDP / MCTS / 归因 / oracle)不感知化学细节。
"""

__all__ = ["__version__"]
__version__ = "0.1.0"
