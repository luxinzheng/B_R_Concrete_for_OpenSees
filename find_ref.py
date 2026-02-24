# -*- coding: utf-8 -*-
import os
root = r'e:\Basic\Concrete_Model\ADINA'
for f in os.listdir(root):
    if f.endswith('.exe') and '\u53c2\u8003' in f:
        full = os.path.join(root, f)
        print(f'Found: {f}')
        print(f'Full path: {full}')
