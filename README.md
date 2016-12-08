# PM-MVS
PM-MVS

The implementation of the PatchMatch MVS in the scene space.

Traditionally, PatchMatch based method needs to traverse all the pixels all the viewpoints several iterations. We instead came up with a way to traverse the 3D points, which significantly improves the speed of the algorithm without losing accuracy.

eval `ssh-agent -s`
ssh-add /c/Users/Admin/.ssh/id_rsa_kw