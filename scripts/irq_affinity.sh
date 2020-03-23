#!/bin/bash
for file in `find /proc/irq -name "smp_affinity"`
do
    var=0x`cat ${file}`
    var="$(( $var & 0x7f ))"
    var=`printf '%.2x' ${var}`
    sudo bash -c "echo ${var} > ${file}"
done
sudo bash -c "echo 7f > /proc/irq/default_smp_affinity"