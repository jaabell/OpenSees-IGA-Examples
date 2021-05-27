#  IGA CANTILEVER PLATE UNDER SELF WEIGHT. ONE ELEMENT MESH, LINEAR CONVERGENCE OBTAINED



set La 10.0   ;# v length
set Lb 1.0    ;# u length
set mm [expr 1.0 / 1000.]  ;# meter to milimeter

wipe
model basic -ndm 3 -ndf  3  ;# 3D model, 3 DOF per node


# These are in u, v order
set controlPts [list \
0 0 0 1 \
0 $Lb 0 1 \
[expr $La*1./4.] 0 0 1 \
[expr $La*1./4.] $Lb 0 1 \
[expr $La*2./4.] 0 0 1 \
[expr $La*2./4.] $Lb 0 1 \
[expr $La*3./4.] 0 0 1 \
[expr $La*3./4.] $Lb 0 1 \
$La 0 0 1 \
$La $Lb 0 1 \
]



# MAterial parameters
set E1 [expr 2.1e11 ] ;# Young's modulus N/m^2
set E2 [expr $E1]
set nu [expr 0 ] ;# Poisson's ratio
set rho [expr 1.0e4]  ;# kg/m^3


# Creating the necessary materials
set tagNDmat1 1
nDMaterial ElasticIsotropic $tagNDmat1 $E1 $nu $rho

set tagNDmat2 2
nDMaterial ElasticIsotropic $tagNDmat2 $E2 $nu $rho


# Creatring the plane stress used for the element formulation
set tagPlaneStress1 3
nDMaterial PlaneStress $tagPlaneStress1 $tagNDmat1

set tagPlaneStress2 4
nDMaterial PlaneStress $tagPlaneStress2 $tagNDmat2

set pi 3.141593
set deg2rad [expr $pi / 180]  ;# for conversion from radian to degrees

set matTags [list $tagPlaneStress1 $tagPlaneStress2 $tagPlaneStress1 $tagPlaneStress2 $tagPlaneStress1]  ;# material tag for each layer
set thickness [list \
 [expr 10. * $mm]\
 [expr 10. * $mm]\
 [expr 10. * $mm]\
 [expr 10. * $mm]\
 [expr 10. * $mm]]  ;# Thickness of each layer
set theta [list \
[expr 0 * $deg2rad]\
[expr 45 * $deg2rad]\
[expr 90 * $deg2rad]\
[expr -45 * $deg2rad]\
[expr 0 * $deg2rad]]  ;# Angle of orientation of each layer

# body force acceleration factors (optional, using selfWeight here)
set gFact [list 0.0  0.0  0.0]

# knot vectors along each direction
set uKnot [list 0 0 1 1]
set vKnot [list 0 0 0 0 0 1 1 1 1 1]


# NURBS basis functions order (P in "u" direction, Q in "v" direction)
set P 1
set Q 4


# Number of control points in each direction
set noPtsX 2
set noPtsY 5

# Number of layers/plies for each IGA patch
set Nlayers [llength theta]


# The patch is an element, so the element tags start with tag "patchTag + 1". Intended for multiPatches
set patchTag  1

# The tag of the first node/control point in the patch (in case of multiPatches, have to update this with the last added node)
set nodeStartTag  1


# ops.IGA call to create the patch (super element)
IGA Patch $patchTag $nodeStartTag $P $Q $noPtsX $noPtsY \
        -type "KLShell" \
        -nonLinearGeometry 1 \
        -planeStressMatTags {*}$matTags \
        -gFact {*}$gFact \
        -theta {*}$theta \
        -thickness {*}$thickness \
        -uKnot {*}$uKnot \
        -vKnot {*}$vKnot \
        -controlPts {*}$controlPts


# Boundary conditions, fixing first two rows for "clamping" condition, fixing "y" displacement for the rest
# for n in ops.getNodeTags():
#     if n in [1, 2, 3, 4]:
#         ops.fix(n, 1, 1, 1)
#     else:
#         ops.fix(n, 0, 1, 0)

foreach node [getNodeTags] {
    if {[lsearch -exact {1 2 3 4} $node] >= 0} {
        puts "Fix $node"
        fix $node 1 1 1
    } else {
        puts "Free $node"
        fix $node 0 1 0
    }
}



puts "\n\n\nPRINTING DOMAIN-----------------------"
printModel
puts "\n\n\nDONE PRINTING DOMAIN-----------------------"


# STKO Recorder

recorder mpco "iga_cantilever" -N displacement -E stresses -E strains 

puts "DONE! "


# exit(0)

# ------------------------------
# Start of analysis generation
# ------------------------------

# create TimeSeries
timeSeries "Linear" 1

# create a plain load pattern
set weight {0.0 0.0 -9.8066}
pattern Plain 1 1 {
    eleLoad -ele  1 -type -SelfWeight  {*}$weight
}

# Loading the patch, selfweight in this case
# -ele 1 means that the patch with tag 1 will load it's associated elements


puts "Starting analysis"

# create SOE
system "FullGeneral"

# create DOF number
numberer "Plain"

# create constraint handler
constraints "Plain"

# create integrator
integrator "LoadControl" 1.0

# create algorithm
test "NormDispIncr" 1.0e-8 50 2
algorithm "Newton"

# create analysis object
analysis "Static"


# perform the analysis
analyze 1

puts "Finished analysis!\n"

puts "Measured vertical displacements at the tip : "
puts "ops.nodeDisp(7,2):  [expr 1000*[nodeDisp 9 3]] mm"
puts "ops.nodeDisp(8,2):  [expr 1000*[nodeDisp 10 3]] mm\n"


proc ladd {l} {::tcl::mathop::+ {*}$l}

set sumThickness [ladd $thickness]

set I [expr ($Lb*pow($sumThickness,3))/12.0]
set W [expr $rho*[lindex $weight 2]*($sumThickness*$Lb)]
set elasticSolution [expr $W*pow($La,4)/(8*$E1*$I)]

puts "elasticSolution:  [expr 1000*$elasticSolution] mm\n"

puts "Done!"

remove "recorders"

exit 0
