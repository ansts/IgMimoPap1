#
# Reassign classes of duplicate peptides that belong to two classes
# Rules applied in this order : 
# -if one copy is from peppos, pep5pred or pepturnd, the others are changed to it too; 
# -if the two copies belong to pepneglo and pepneg both should be pepneglo;
# -if the two copies belong to pepother5 and pepneg both should be pepother5
#
reassdup=function(duppp){
dupreass=list()
dupp=unique(duppp$pep)
for (x in dupp){
  cxx=duppp$pep==x
  if ("peppos" %in% duppp$class[cxx]) {
    frag=duppp[cxx,]
    frag[frag[,1]!="peppos",1]="peppos"
    dupreass=rbind(dupreass,frag)
    next}
  if ("pep5pred" %in% duppp$class[cxx]) {
    frag=duppp[cxx,]
    frag[frag[,1]!="pep5pred",1]="pep5pred"
    dupreass=rbind(dupreass,frag)
    next}
  if ("pepturnd" %in% duppp$class[cxx]) {
    frag=duppp[cxx,]
    frag[frag[,1]!="pepturnd",1]="pepturnd"
    dupreass=rbind(dupreass,frag)
    next}
  if (("pepneglo" %in% duppp$class[cxx]) & ("pepneg" %in% duppp$class[cxx])) {
    frag=duppp[cxx,]
    frag[frag[,1]!="pepneglo",1]="pepneglo"
    dupreass=rbind(dupreass,frag)
    next}
  if (("pepother5" %in% duppp$class[cxx]) & ("pepneg" %in% duppp$class[cxx])) {
    frag=duppp[cxx,]
    frag[frag[,1]!="pepother5",1]="pepother5"
    dupreass=rbind(dupreass,frag)
    next}
}
return(dupreass)
}