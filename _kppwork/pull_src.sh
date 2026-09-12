#!/bin/bash
D=/mnt/d/atmos/_kppwork/src
cp /root/kpp2/cb05_Parameters.f90 $D/old_Parameters.f90
cp /root/kpp2/cb05_Monitor.f90    $D/old_Monitor.f90
cp /root/kpp2/cb05_Function.f90   $D/old_Function.f90
cp /root/kpp5/cb05_Parameters.f90 $D/new_Parameters.f90
cp /root/kpp5/cb05_Monitor.f90    $D/new_Monitor.f90
cp /root/kpp5/cb05_Function.f90   $D/new_Function.f90
ls $D | grep -E "^(old|new)_"
