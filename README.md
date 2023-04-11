# SoilLab
Tools for soil physics lab

## Identification file
Contains metadata for all samples in a project: The lab ID, ring number, sample metadata, analysis and comments.
It is a tab delimited text file.
Number of headerlines are standard = 6
The data has 2 column headers: the first one is used to seperate cathegories of data, the second one specifies them. These headers should remain mostly untouched, but the 'Metadata' section is flexible. Columns can be added or removed, empty columns are ignored. The MetaLabels can be modified. The last metalabel is assumed to be the sample number.

Example:
```
Project	Test													
Date	10/05/17
Collaborators	Jan	 Maarten

Identification			Metadata					Analysis						Comments
LabID	BatchID	RingID	MetaLabel1	MetaLabel2	MetaLabel3	MetaLabel4	MetaLabel5	SBPP	Permeameter	AirPermeability	UMS_Hyprop	UMS_Ksat	Texture	Comments
17.054	1	654	Field1	Top	Head		654					1		
17.054	1	321	Field1	Sub	Head		321					1		
17.054	1	984	Field1	Top	Center		984					1		
17.054	1	987	Field1	Sub	Center		987					1		
```

## SBPP (Sandbox Pressure Plates)

Example:
```
Project	Test																								
Date	10/05/17																								
Collaborators	Jan, Maarten																								
																									
Identification			Metadata			Comments	Data																		
LabID	BatchID	RingID	MetaLabel1	MetaLabel3	MetaLabel4	Comments	RingHeight(m)	RingDiameter(m)	EmptyRingWeight(kg)	NylonWeight(kg)	MassStable10(kg)	MassStable50(kg)	MassStable100(kg)	MassSubsample100Cup(kg)	MassSubsample100Wet(kg)	MassSubsample100Dry(kg)	MassSubsample340Cup(kg)	MassSubsample340Wet(kg)	MassSubsample340Dry(kg)	MassSubsample1020Cup(kg)	MassSubsample1020Wet(kg)	MassSubsample1020Dry(kg)	MassSubsample15400Cup(kg)	MassSubsample15400Wet(kg)	MassSubsample15400Dry(kg)
17.013	1	15	Field1	Head	15		0.049847	0.050127	0.099846	0.0019984	0.23946	0.22608	0.21714	0.024057	0.037804	0.035844	0.011532	0.031347	0.029555	0.011508	0.031381	0.029379	0.011475	0.031304	0.029607
17.013	1	16	Field1	Head	16		0.050168	0.050198	0.10008	0.0020017	0.24114	0.22431	0.21771	0.023978	0.037514	0.036089	0.01144	0.031236	0.029511	0.011452	0.031223	0.029451	0.011469	0.031232	0.029512
```

## UMS_Hyprop

All exported data in a folder, given a proper name. Calculations and Fit need to be executed in Hyprop.

## UMS_Ksat

All data in a folder, given a proper name.
