FROM rocker/shiny

RUN R -e "install.packages(c('ggplot2','annotate','Biobase','data.table','EMA','GEOmetadb','GEOquery','ggplot2','gplots','knitr','latex2exp','limma','mygene','R.utils'), repos='http://cran.rstudio.com/')"

COPY R/test_App /srv/shiny-server/test_app

CMD ["/usr/bin/shiny-server.sh"]