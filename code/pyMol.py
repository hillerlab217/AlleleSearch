reinit
### Retrieving the PDB file ###
fetch 1395, Kchannel, async=0
### Establishing selections and objects ###
create SFilter, resi 75-80
select K, name K
deselect
###/Visualisation ###
remove resn TBA
util.cbc
as sticks, SFilter color white, SFilter util.cbaw SFilter util.cbaw K
###Orientation ###
###Savescene ###
scene F1, store
###New scene ###
hide sticks, SFilter
create Kcomplex, resi 77+78+201
as sticks, Kcomplex
util.cbaw Kcomplex
hide sphere, resi 202-204