from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors

#Writes map to a JPG
def WriteImg(seqRecIn):
    #to select different colors from
    colorSelection = [colors.darkgreen, colors.darkorange, colors.lightblue]
    currColor = 0
    #uses GenomeDiagram package to create a featureset object for the diagram to use later
    featSet = GenomeDiagram.FeatureSet()
    for currFeat in seqRecIn.features:
        if currFeat.type == "CDS":
            featSet.add_feature(currFeat,
            color = colorSelection[currColor%3],
            label = True,
            label_size=20,
            sigil="ARROW",
            #alternates widths specifically to accomodate the TCSV Sequence: order goes Medium->Thinnest->Thickest->Medium->Thinnest
            arrowshaft_height = (0.5 + ((0.2) - (((currColor+1)%3) * (0.2)))),
            #turns arrowhead triangle into a thin line
            arrowhead_length = 0)
            currColor += 1
    #add features to track(only 1 because it is circular)
    track = GenomeDiagram.Track(name="TCSV Sequence Features")
    track.add_set(featSet)
    #add track to diagram
    diagram = GenomeDiagram.Diagram(name="Tomato Curly Stunt Virus",
    format="circular",
    pagesize=(1000,1000),
    circular=True,
    start=0,
    #goes until end of any length sequence
    end=len(seqRecIn),
    circle_core=0.5)
    diagram.add_track(track,1)
    #writes diagram to memory so it can be written into a file
    diagram.draw()
    #writes to a file using the reportlab functionalities
    diagram.write("TCSV_Sequence_Map.jpg", "JPG")

#parses file
def ReadSeq(fileName, fileFormat):
    #assumes there will only be one sequence in the file
    seqRec = SeqIO.read(fileName, fileFormat)
    WriteImg(seqRec)

ReadSeq("Genome.gb", "genbank")