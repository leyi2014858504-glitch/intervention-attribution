#!/bin/bash
python3 -c "import sklearn; print('wsl sklearn', sklearn.__version__)" 2>&1 | tail -1