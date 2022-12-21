# ShinyRFU
The ShinyRFU() program is in file: ShinyRFU_0.2.0.tar.gz (latest version uploaded 3/11/2022). The program itself is written in R and it utilises EuroForMix http://www.euroformix.com/ in the background in order to carry out the calculations. There is a user friendly gui interface incorporated, that is programmed in R-Shiny https://shiny.rstudio.com/. It is not intended as a replacement for EuroForMix as it has no exploratory capability: it is primarily intended to provide a quick method to calculate average RFU value and associated mixture proportion (Mx) values. If there are a large number of samples to analyse, they can be concatenated into a single file  and analysed in a single run – this is a considerable time saving. 

Use v.0.2.0. Previous version v.0.1.1 has a bug that caused the program to crash when conditioning information is absent

Tested with R4.2.2
Make sure you load latest version of EUroformix 4.0.1
