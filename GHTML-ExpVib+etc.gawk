#!/usr/bin/awk -f
# parse Expt@cccbdb HTML 
# stdin should be handled specially, due to FILENAME processing!
#BEWARE Win version does not seem to work?
# script needs to have system-specific CR formatting
#BEWARE use HTML, as raw output table does not distinguish columns _
#	Fundamental/Harmonic in frequency table!

# this hangs with input filename that contains "="!

BEGIN {
#print "TEST delSubSup(text)='" delSubSup("<td><sup>1</sup>A<sub>1g</sub></td>") "'"
	OFS="\t"
	nProp=0
	tProp[++nProp]="inchi"	# tags (string labels/indexes) for properties
	tProp[++nProp]="SMILES"	# tags (string labels/indexes) for properties
# not needed:	iProp["inchi"]=nProp; numindex for tags (string labels/indexes) for properties
	tProp[++nProp]="Electronic State"
	tProp[++nProp]="'Conformation'"
	tProp[++nProp]="Hfg(298.15K)_val"
	tProp[++nProp]="Hfg(298.15K)_unit"
	tProp[++nProp]="Hfg(0K)_val"
	tProp[++nProp]="Hfg(0K)_unit"
	tProp[++nProp]="Entropy(298.15K)_val"
	tProp[++nProp]="Entropy(298.15K)_unit"
	tProp[++nProp]="IntegratedHeatCapacity(0to298.15K)_val"
	tProp[++nProp]="IntegratedHeatCapacity(0to298.15K)_unit"
	tProp[++nProp]="HeatCapacity(298.15K)_val"
	tProp[++nProp]="HeatCapacity(298.15K)_unit"
	tProp[++nProp]="ZPVE"
	tProp[++nProp]="dVib1Line"
	tProp[++nProp]="NIST_ID"
	tProp[++nProp]="Point_Group"
	 
# property values are to be collected in vProp[tag] array
# currect code only works with single-file input, so this vProp[tag] is 1D
#print "TEST extrTagged(text, num): " extrTagged("<1><2>after 2<3>", 2)

	if (!selectVF) selectVF="Fundamental"
	#selectVF:
	#tVibHead[2,1]="Fundamental"
	#tVibHead[2,2]="Harmonic"
	if (selectVF=="Fundamental") {
		iVibCol=3
	} else {
		iVibCol=4
	}
}

(FNR==1 && NR>1) {
	print "ERROR: this script only works on single file input"; exit 
}

(FNR==1) {
	filebase=stripext(FILENAME)
#	print "NR="NR, filebase
	section=""
	subsection=""
}

/<table /{
	table=1
}

/<\/table /{
	table=0
	section=""
}

(table)&&($0 ~ /<th>State/){
	section="State"
	iTR=0
}

(table)&&(section=="State")&&($0 ~ /<th>Conformation/){
	subsection="State&Conformation"
	iTD=0
}

(table)&&(subsection=="State&Conformation")&&($0 ~ /<td/){
	iTD++
	if (iTD==1) {
		#print "track Electronic State: "$0
		vProp["Electronic State"]=extrTagged(delSubSup($0), 1)
	}
	if (iTD==2) {
		vProp["'Conformation'"]=extrTagged($0, 1)
		subsection=""
#		print "Electronic State", vProp["Electronic State"]
#		print "'Conformation'", vProp["'Conformation'"]
	}
}

/<th.*International Chemical Identifier/{
	section="inchi+etc"
	iTD=0
}

(section=="inchi+etc")&&($0 ~ /<td>/){
	iTD++
	if (iTD==1) {
		nextData="inchi"
		vProp[nextData]=extrTagged($0, 2)
#		print nextData, vProp[nextData]
	}
	if (iTD==3) {
		nextData="SMILES"
		vProp[nextData]=extrTagged($0, 1)
#		print nextData, vProp[nextData]
	}
	if (iTD==4) {
		nextData="IUPAC"
		vProp[nextData]=extrTagged($0, 2)
#		print nextData, vProp[nextData]
		section=""
	}
}


($0 ~ /<div.*title="Enthalpy and Entropy/){
	section="Enthalpy+etc"
	iTR=0
}

(table)&&(section=="Enthalpy+etc")&&($0 ~ /<th>Property/){
	subsection="Enthalpy+etc_Data"
	iTR=0
}

(section == "Enthalpy+etc")&&(subsection == "Enthalpy+etc_Data")&&($0 ~ /<tr/){
	#print "TRACE: @ Enthalpy+etc @: " $0
	iTR++
	iTD=0
}

(section == "Enthalpy+etc")&&(subsection == "Enthalpy+etc_Data")&&($0 ~ /<td/){
	#print "TRACE: @ Enthalpy+etc @: " $0
	iTD++
#print "TRACK " subsection " - iTD="iTD", $0: " $0
	if (iTD==1) {
		nextLabel=extrTagged($0, 1) 
		gsub(" *", "", nextLabel)
	}
	if (iTD==2) {
		vProp[nextLabel "_val"]=extrTagged($0, 1) 
	}
	if (iTD==4) {
		vProp[nextLabel "_unit"]=extrTagged(delSubSup($0), 1)
#		print nextLabel, vProp[nextLabel "_val"], vProp[nextLabel "_unit"] 
	}
	
}

###

/Vibrational symmetries, frequencies, and intensities/{
	section="frequencies"
	subsection=""
	vProp["dVib1Line"]=""
}

/vibrational zero-point energy:/{ # cm<SUP>-1</SUP> (from fundamental vibrations)
	nDVib=iTR # this closed the "frequencies" section, so assign their data count
#	print "dVib1Line", dVib1Line 
	section="zero-point energy"
	subsection=""
	vProp["ZPVE"]=$0
	vProp["ZPVE"]=floatNumberExtract(vProp["ZPVE"])
#	print "ZPVE", vProp["ZPVE"]
	section=""
}

(section == "frequencies")&&($0 ~ /<table /){
	#print "TRACE: @ frequencies @: " $0
	subsection="header"
	iTR=0
}

(section == "frequencies")&&(subsection == "header")&&($0 ~ /<tr/){
	#print "TRACE: @ frequencies @: " $0
	iTR++
	iTH=0
	if (iTR>2) {
		subsection="vibData"
		iTR=0
	}
}

(section == "frequencies")&&(subsection == "header")&&($0 ~ /<th/){
	#print "TRACE: @ frequencies @: " $0
	iTH++
	tVibHead[iTR,iTH]=extrTagged($0, 1)
	#print "TRACK header@frequencies: tVibHead["iTR","iTH"]=" tVibHead[iTR,iTH]
	#selectVF:
	#tVibHead[2,1]="Fundamental"
	#tVibHead[2,2]="Harmonic"
}

(section == "frequencies")&&(subsection == "vibData")&&($0 ~ /<tr/){
	#print "TRACE: @ frequencies @: " $0
	iTR++
	iTD=0
}

(section == "frequencies")&&(subsection == "vibData")&&($0 ~ /<td/){
	nextLabel="dVib1Line"
	iTD++
	if (iTD==1) {
		dVibData[iTR,1]=extrTagged($0, 1) 
		if (iTR != dVibData[iTR,1]) print FILENAME " WARNING data line mismatch at iTR="iTR", dVibData[iTR,1]="dVibData[iTR,1]
	}
	if (iTD==2) {
		sub("<sub>","[")
		sub("^[^>]*>","")
		sub("<.*$","]")
		dVibData[iTR,2]=extrTagged(delSubSup($0), 1) 
	}
#selectVF:iVibCol
	if (iTD==iVibCol) {
		dVibData[iTR,3]=extrTagged($0, 1) 
		vProp[nextLabel]=vProp[nextLabel] OFS dVibData[iTR,3]
	}
}
(section == "frequencies")&&(subsection == "vibData")&&($0 ~ /<\/table/){
	nVibData=iTR
	iTR=0
	section=""
	subsection=""
}

###
/cbook\.cgi\?ID=/{
	nextLabel="NIST_ID"
	sub("^.*cbook.cgi\?ID=", "")
	sub("&amp;.*$", "")
	sub("\".*$", "")
	vProp[nextLabel]=$0
#	print nextLabel, vProp[nextLabel] 
	nextLabel=""
}

/<P.*Point Group /{
	nextLabel="Point_Group"
#<P class="strong">Point Group D<sub>6h</sub></p>
	sub("Point Group ", "")
	vProp[nextLabel]=extrTagged(delSubSup($0), 1)
	#print "TRACK Point_Group: nextLabel='"nextLabel"'", vProp[nextLabel] 
#	print nextLabel, vProp[nextLabel] 
	nextLabel=""
}

($0 ~ /<br.*Internal coordinates/){
	section="IntCoo"
	iTR=0
}

(table)&&(section=="IntCoo")&&($0 ~ /<th/){
	subsection="IntCoo_Head"
	iTR=0
}


(section=="IntCoo")&&(subsection=="IntCoo_Head")&&($0 ~ /<td/){
	subsection="IntCoo_Tab"
	iTR=0
}

(section == "IntCoo")&&(subsection == "IntCoo_Tab")&&($0 ~ /<tr/){
	iTR++
	iTD=0
}

(section == "IntCoo")&&(subsection == "IntCoo_Tab")&&($0 ~ /<td/){
	#print "TRACE: @ Enthalpy+etc @: " $0
	iTD++
#print "TRACK " subsection " - iTD="iTD", $0: " $0
	if (1<=iTD&& iTD<=2) {
		aIntCoo[iTR,iTD]=extrTagged($0, 1) 
	}
	if (4<=iTD&& iTD<=7) {
		aIntCoo[iTR,iTD]=extrTagged($0, 1) 
	}
}
(section == "IntCoo")&&(subsection == "IntCoo_Tab")&&($0 ~ /<\/table/){
	nIntCoo=iTR
	iTR=0
	section=""
	subsection=""
}

/<span.*Cartesians/{
	section="CartCoo"
	iTR=0
}

(table)&&(section=="CartCoo")&&($0 ~ /<th/){
	subsection="CartCoo_Tab"
	iTR=0
}

(section == "CartCoo")&&(subsection == "CartCoo_Tab")&&($0 ~ /<tr/){
	iTR++
	iTD=0
}

(section == "CartCoo")&&(subsection == "CartCoo_Tab")&&($0 ~ /<td/){
	iTD++
#print "TRACK " subsection " - iTD="iTD", $0: " $0
	if (iTD==1) {
		aCartCoo[iTR,iTD]=extrTagged($0, 2) 
	}
	if (2<=iTD&& iTD<=4) {
		aCartCoo[iTR,iTD]=extrTagged($0, 1) 
	}
}
(section == "CartCoo")&&(subsection == "CartCoo_Tab")&&($0 ~ /<\/table/){
	nCartCoo=iTR
	iTR=0
	section=""
	subsection=""
	# end of processing, print output files
	output(filebase)
}

function stripext(file) {
	sub("\.[^\.]*$", "", file)
	return file
}

function floatNumberExtract(text) {# extract the 1st number within text
	text=substr(text, match(text, "[0-9][0-9.]+"), RLENGTH)
#	text=substr(text, 1, index(text, "[^0-9.]"))
	return text
}

function extrTagged(text, num) {# extract the field after 'num' HTML tags 
	for (i=1; i<=num; i++) {
#print "DEBUG extrTagged: "text
#print "DEBUG extrTagged: RLENGTH="RLENGTH
		match(text, "^[^>]*>")
		text=substr(text, RLENGTH+1)
	}
	sub("<.*$", "", text)
#	text=substr(text, 1, index(text, "[^0-9.]"))
	return text
}

#<td><sup>1</sup>A<sub>1g</sub></td>
function delSubSup(text) {
	gsub("<sup>", "", text)
	gsub("</sup>", "", text)
	gsub("<sub>", "", text)
	gsub("</sub>", "", text)
#print "DEBUG delSubSup: text='"text"'"
	return text
}

function output(filebase) {
	propFile = filebase ".prop"
	xyzFile = filebase ".xyz"
	mopFile = filebase ".mop"
	vibFile = filebase "_Fundamentals.csv"
	for (iProp=1; iProp<=nProp; iProp++) print tProp[iProp], vProp[tProp[iProp]] >propFile

	print nCartCoo >xyzFile
	print vProp["SMILES"] vProp["dVib1Line"] >xyzFile
	for (i=1; i<=nCartCoo; i++) {
		out = aCartCoo[i,1]
		gsub("[0-9]", "", out)
		for (j=2; j<=4; j++) {
			out = out " " aCartCoo[i,j]
		}
		print out>xyzFile
	}
	print "">xyzFile

	print "PM6 SYMMETRY" >mopFile
	print vProp["SMILES"] vProp["dVib1Line"] >mopFile
	print vProp["inchi"] >mopFile
	for (i=1; i<=nCartCoo; i++) {
		out = aCartCoo[i,1]
		gsub("[0-9]", "", out)
		for (j=2; j<=4; j++) {
			out = out " " aCartCoo[i,j] " +1"
		}
		print out>mopFile
	}
	print " ">mopFile #BEWARE obabel chokes on empty blank line!

	for (i=1; i<=nVibData; i++) {
		out = dVibData[i,1]
		for (j=2; j<=3; j++) {
			out = out "," dVibData[i,j] 
		}
		print out>vibFile
	}
	print "">vibFile

#dVibData[iTR,3]
}

