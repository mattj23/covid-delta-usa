from typing import Dict, List, Tuple


def state_to_abbrev() -> Dict[str, str]:
    return {n: a for n, a in _split_raw()}


def abbrev_to_state() -> Dict[str, str]:
    return {a: n for n, a in _split_raw()}


def state_abbrevs() -> List[str]:
    return [a for _, a in _split_raw()]


def _split_raw() -> List[Tuple[str, str]]:
    pairs = []
    for line in _raw_data.split("\n"):
        if not line.strip():
            continue

        name, abbrev, *_ = line.split("\t")
        pairs.append((name, abbrev))
    return pairs


def adjacency(state: str) -> List[str]:
    return _compute_adjacencies()[state]


def _compute_adjacencies() -> Dict[str, List[str]]:
    clusters = [x.split(",") for x in _adjacency.split("\n") if x.strip()]
    results = {}
    for c in clusters:
        results[c[0]] = c[1:]
    return results


_adjacency = """AK,WA,OR
AL,MS,TN,GA,FL
AR,MO,TN,MS,LA,TX,OK
AZ,CA,NV,UT,CO,NM
CA,OR,NV,AZ
CO,WY,NE,KS,OK,NM,AZ,UT
CT,NY,MA,RI
DE,MD,PA,NJ
FL,AL,GA
GA,FL,AL,TN,NC,SC
HI,CA
IA,MN,WI,IL,MO,NE,SD
ID,MT,WY,UT,NV,OR,WA
IL,IN,KY,MO,IA,WI
IN,MI,OH,KY,IL
KS,NE,MO,OK,CO
KY,IN,OH,WV,VA,TN,MO,IL
LA,TX,AR,MS
MA,RI,CT,NY,NH,VT
MD,VA,WV,PA,DC,DE
ME,NH
MI,WI,IN,OH
MN,WI,IA,SD,ND
MO,IA,IL,KY,TN,AR,OK,KS,NE
MS,LA,AR,TN,AL
MT,ND,SD,WY,ID
NC,VA,TN,GA,SC
ND,MN,SD,MT
NE,SD,IA,MO,KS,CO,WY
NH,VT,ME,MA
NJ,DE,PA,NY
NM,AZ,UT,CO,OK,TX
NV,ID,UT,AZ,CA,OR
NY,NJ,PA,VT,MA,CT
OH,PA,WV,KY,IN,MI
OK,KS,MO,AR,TX,NM,CO
OR,CA,NV,ID,WA
PA,NY,NJ,DE,MD,WV,OH
RI,CT,MA
SC,GA,NC
SD,ND,MN,IA,NE,WY,MT
TN,KY,VA,NC,GA,AL,MS,AR,MO
TX,NM,OK,AR,LA
UT,ID,WY,CO,NM,AZ,NV
VA,NC,TN,KY,WV,MD,DC
VT,NY,NH,MA
WA,ID,OR
WI,MI,MN,IA,IL
WV,OH,PA,MD,VA,KY
WY,MT,SD,NE,CO,UT,ID
"""

_raw_data = """
Alabama	AL	Ala.
Alaska	AK	Alaska
Arizona	AZ	Ariz.
Arkansas	AR	Ark.
California	CA	Calif.
Colorado	CO	Colo.
Connecticut	CT	Conn.
Delaware	DE	Del.
Florida	FL	Fla.
Georgia	GA	Ga.
Hawaii	HI	Hawaii
Idaho	ID	Idaho
Illinois	IL	Ill.
Indiana	IN	Ind.
Iowa	IA	Iowa
Kansas	KS	Kans.
Kentucky	KY	Ky.
Louisiana	LA	La.
Maine	ME	Maine
Maryland	MD	Md.
Massachusetts	MA	Mass.
Michigan	MI	Mich.
Minnesota	MN	Minn.
Mississippi	MS	Miss.
Missouri	MO	Mo.
Montana	MT	Mont.
Nebraska	NE	Neb. or Nebr.
Nevada	NV	Nev.
New Hampshire	NH	N.H.
New Jersey	NJ	N.J.
New Mexico	NM	N.Mex.
New York	NY	N.Y.
North Carolina	NC	N.C.
North Dakota	ND	N.Dak.
Ohio	OH	Ohio
Oklahoma	OK	Okla.
Oregon	OR	Ore. or Oreg.
Pennsylvania	PA	Pa.
Rhode Island	RI	R.I.
South Carolina	SC	S.C.
South Dakota	SD	S.Dak.
Tennessee	TN	Tenn.
Texas	TX	Tex. or Texas
Utah	UT	Utah
Vermont	VT	Vt.
Virginia	VA	Va.
Washington	WA	Wash.
West Virginia	WV	W.Va.
Wisconsin	WI	Wis. or Wisc.
Wyoming	WY	Wyo.
"""

if __name__ == '__main__':
    _compute_adjacencies()
