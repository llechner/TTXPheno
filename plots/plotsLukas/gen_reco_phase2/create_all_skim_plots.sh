#!/bin/bash

cpQM=7
cpt=8
ctW=6
ctWI=6
ctZ=6
ctZI=6
ctG=6
ctGI=6

#Run all skim plot scripts

#./skim_plots_ttgamma_1l.sh ${cpQM} ${cpt} ${ctW} ${ctWI} ${ctZ} ${ctZI} #${ctG} ${ctGI}
#./skim_plots_ttgamma_2l.sh ${cpQM} ${cpt} ${ctW} ${ctWI} ${ctZ} ${ctZI} ${ctG} ${ctGI}
./skim_plots_ttZ_3l.sh ${cpQM} ${cpt} ${ctW} ${ctWI} ${ctZ} ${ctZI} #${ctG} ${ctGI}
#./skim_plots_ttZ_4l.sh ${cpQM} ${cpt} ${ctW} ${ctWI} ${ctZ} ${ctZI} #${ctG} ${ctGI}
#./skim_plots_ttW_2l.sh ${cpQM} ${cpt} ${ctW} ${ctWI} ${ctZ} ${ctZI} #${ctG} ${ctGI}

exit

cpQM=5
cpt=5
ctW=5
ctWI=5
ctZ=5
ctZI=5
ctG=3
ctGI=3

#Run all skim plot scripts

./skim_plots_ttZ_3l.sh ${cpQM} ${cpt} ${ctW} ${ctWI} ${ctZ} ${ctZI} ${ctG} ${ctGI}
./skim_plots_ttZ_4l.sh ${cpQM} ${cpt} ${ctW} ${ctWI} ${ctZ} ${ctZI} ${ctG} ${ctGI}
./skim_plots_ttgamma_1l.sh ${cpQM} ${cpt} ${ctW} ${ctWI} ${ctZ} ${ctZI} ${ctG} ${ctGI}
./skim_plots_ttgamma_2l.sh ${cpQM} ${cpt} ${ctW} ${ctWI} ${ctZ} ${ctZI} ${ctG} ${ctGI}
#./skim_plots_ttW_2l.sh ${cpQM} ${cpt} ${ctW} ${ctWI} ${ctZ} ${ctZI} ${ctG} ${ctGI}



cpQM=10
cpt=10
ctW=5
ctWI=5
ctZ=10
ctZI=10
ctG=3
ctGI=3

#Run all skim plot scripts

./skim_plots_ttZ_3l.sh ${cpQM} ${cpt} ${ctW} ${ctWI} ${ctZ} ${ctZI} ${ctG} ${ctGI}
./skim_plots_ttZ_4l.sh ${cpQM} ${cpt} ${ctW} ${ctWI} ${ctZ} ${ctZI} ${ctG} ${ctGI}
./skim_plots_ttgamma_1l.sh ${cpQM} ${cpt} ${ctW} ${ctWI} ${ctZ} ${ctZI} ${ctG} ${ctGI}
./skim_plots_ttgamma_2l.sh ${cpQM} ${cpt} ${ctW} ${ctWI} ${ctZ} ${ctZI} ${ctG} ${ctGI}
#./skim_plots_ttW_2l.sh ${cpQM} ${cpt} ${ctW} ${ctWI} ${ctZ} ${ctZI} ${ctG} ${ctGI}
