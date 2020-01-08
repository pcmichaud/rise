
* directory for datasource
global src "data"
* make a temporary directory to save temp datasets
! mkdir -p $src/temp

* global graph options
global fig "graphregion(color(white)) bgcolor(white)"
* display time stamp
display "$S_TIME  $S_DATE"
* install requried packages (dm83)
set processors 4
