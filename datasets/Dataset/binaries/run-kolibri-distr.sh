#!/bin/sh
qemu-system-x86_64 \
    -name 'Kolibrios' \
    -m 10000\
    -enable-kvm \
    -snapshot \
    -net nic,model=virtio \
    -net user,hostfwd=tcp:127.0.0.1:1030-:22 \
    -serial mon:stdio \
    -fda kolibri.img -boot a
