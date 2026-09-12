"""化学层:可替换的箱模型后端。

- `base`  : `ChemistryBackend` 抽象契约(状态 = 后端全物种)
- `mock`  : 教学用测试替身(编码 VOC-limited 非单调性,毫秒级)
- `cb05`  : 真实气相机制后端(经 pykpp 解析 CB05,79 物种)
"""
from .base import ChemistryBackend, StepResult
from .mock import MockBackend
from .cb05 import CB05Backend

__all__ = ["ChemistryBackend", "StepResult", "MockBackend", "CB05Backend"]
