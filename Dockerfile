USER root

RUN apt-get update && \
	apt-get install -y git
 
ENV HOME /home/developer
WORKDIR $HOME

USER developer

RUN pip install --user --no-cache-dir notebook==5.*
RUN pip install --user future
RUN pip install --user pandas
RUN pip install --user scipy
RUN pip install --user numpy
RUN pip install --user matplotlib
RUN pip install --user ipykernel==4.7.0
RUN pip install --user pygfunction

COPY Introduction.ipynb $HOME
COPY PartI_TABS.ipynb $HOME
COPY PartII_StaticCalc.ipynb $HOME
COPY PartIII_EnergySystemDesign.ipynb $HOME
COPY PartIV_DynamicSim.ipynb $HOME
COPY PartV_hybridGEOTABS.ipynb $HOME
COPY loadCalc.xlsx $HOME
COPY resistanceCalculator.py $HOME
COPY fig $HOME/fig
COPY results $HOME/results
