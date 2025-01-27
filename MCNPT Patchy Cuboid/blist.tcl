set blist {}
set natoms [molinfo top get numatoms]
#for {set a 1; set b 2} {$b < $natoms} {incr a; incr b} {
  # each bond list entry is: <idx> <idx> [<type>] [<order>]
#  lappend blist [list $a $b A-A 1]
#  lappend blist [list $a $b A-A 1]
#}
topo clearbonds
topo getbondlist
#set b "19" 
# for full cube
set b "19+1"
#for vertices and center alone
for {set a 0} {$a < 500} {incr a} {
	set i [expr $a*$b+1]; set j [expr $a*$b+2];
	lappend blist [list $i $j A-A 1]
#	puts $i
#	puts $j
	set i [expr $a*$b+1]; set j [expr $a*$b+3];
	lappend blist [list $i $j A-A 1]
	set i [expr $a*$b+1]; set j [expr $a*$b+5];
	lappend blist [list $i $j A-A 1]
	set i [expr $a*$b+2]; set j [expr $a*$b+4];
	lappend blist [list $i $j A-A 1]
	set i [expr $a*$b+2]; set j [expr $a*$b+6];
	lappend blist [list $i $j A-A 1]
	set i [expr $a*$b+3]; set j [expr $a*$b+4];
	lappend blist [list $i $j A-A 1]
	set i [expr $a*$b+3]; set j [expr $a*$b+7];
	lappend blist [list $i $j A-A 1]
	set i [expr $a*$b+4]; set j [expr $a*$b+8];
	lappend blist [list $i $j A-A 1]
	set i [expr $a*$b+5]; set j [expr $a*$b+6];
	lappend blist [list $i $j A-A 1]
	set i [expr $a*$b+5]; set j [expr $a*$b+7];
	lappend blist [list $i $j A-A 1]
	set i [expr $a*$b+6]; set j [expr $a*$b+8];
	lappend blist [list $i $j A-A 1]
	set i [expr $a*$b+7]; set j [expr $a*$b+8];
	lappend blist [list $i $j A-A 1]

}
topo setbondlist both $blist
vmdcon -info "assigned [topo numbondtypes] bond types to [topo numbonds] bonds:"
vmdcon -info "bondtypes: [topo bondtypenames]"
