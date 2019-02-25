#!/bin/bash

python weight_function_plots.py --small --sample fwlite_ttZ_ll_LO_order2_15weights --parameters ctp ctpI cpQM cpQ3 cpt cptb cptbI ctW ctWI ctZ ctZI cbW cbWI ctG ctGI
python weight_function_plots.py --small --sample fwlite_ttZ_ll_LO_order2_15weights --parameters ctp ctpI cpQM cpQ3 cpt cptb cptbI ctW ctWI ctZ ctZI cbW cbWI ctG ctGI --differential
python weight_function_plots.py --small --sample fwlite_ttZ_ll_LO_order2_15weights --parameters ctp ctpI cpQM cpQ3 cpt cptb cptbI ctW ctWI ctZ ctZI cbW cbWI ctG ctGI --normalized

python weight_function_plots.py --small --sample fwlite_ttZ_ll_LO_order2_15weights_ref --parameters ctp ctpI cpQM cpQ3 cpt cptb cptbI ctW ctWI ctZ ctZI cbW cbWI ctG ctGI
python weight_function_plots.py --small --sample fwlite_ttZ_ll_LO_order2_15weights_ref --parameters ctp ctpI cpQM cpQ3 cpt cptb cptbI ctW ctWI ctZ ctZI cbW cbWI ctG ctGI --differential
python weight_function_plots.py --small --sample fwlite_ttZ_ll_LO_order2_15weights_ref --parameters ctp ctpI cpQM cpQ3 cpt cptb cptbI ctW ctWI ctZ ctZI cbW cbWI ctG ctGI --normalized

python weight_function_plots.py --small --sample fwlite_ttZ_ll_LO_order3_8weights --order 3 --parameters cpQM cpt ctW ctWI ctZ ctZI ctG ctGI
python weight_function_plots.py --small --sample fwlite_ttZ_ll_LO_order3_8weights --order 3 --parameters cpQM cpt ctW ctWI ctZ ctZI ctG ctGI --differential
python weight_function_plots.py --small --sample fwlite_ttZ_ll_LO_order3_8weights --order 3 --parameters cpQM cpt ctW ctWI ctZ ctZI ctG ctGI --normalized

python weight_function_plots.py --sample fwlite_ttZ_ll_LO_order2_15weights --parameters ctp ctpI cpQM cpQ3 cpt cptb cptbI ctW ctWI ctZ ctZI cbW cbWI ctG ctGI
python weight_function_plots.py --sample fwlite_ttZ_ll_LO_order2_15weights --parameters ctp ctpI cpQM cpQ3 cpt cptb cptbI ctW ctWI ctZ ctZI cbW cbWI ctG ctGI --differential
python weight_function_plots.py --sample fwlite_ttZ_ll_LO_order2_15weights --parameters ctp ctpI cpQM cpQ3 cpt cptb cptbI ctW ctWI ctZ ctZI cbW cbWI ctG ctGI --normalized

python weight_function_plots.py --sample fwlite_ttZ_ll_LO_order2_15weights_ref --parameters ctp ctpI cpQM cpQ3 cpt cptb cptbI ctW ctWI ctZ ctZI cbW cbWI ctG ctGI
python weight_function_plots.py --sample fwlite_ttZ_ll_LO_order2_15weights_ref --parameters ctp ctpI cpQM cpQ3 cpt cptb cptbI ctW ctWI ctZ ctZI cbW cbWI ctG ctGI --differential
python weight_function_plots.py --sample fwlite_ttZ_ll_LO_order2_15weights_ref --parameters ctp ctpI cpQM cpQ3 cpt cptb cptbI ctW ctWI ctZ ctZI cbW cbWI ctG ctGI --normalized

python weight_function_plots.py --sample fwlite_ttZ_ll_LO_order3_8weights --order 3 --parameters cpQM cpt ctW ctWI ctZ ctZI ctG ctGI
python weight_function_plots.py --sample fwlite_ttZ_ll_LO_order3_8weights --order 3 --parameters cpQM cpt ctW ctWI ctZ ctZI ctG ctGI --differential
python weight_function_plots.py --sample fwlite_ttZ_ll_LO_order3_8weights --order 3 --parameters cpQM cpt ctW ctWI ctZ ctZI ctG ctGI --normalized
