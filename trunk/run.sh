#!/bin/sh
echo bfilter
#nice -19 bin/bfilter $1 $2 $3 $4 
echo bfilterOMP
#nice -19 bin/bfilterOMP $1 $2 $3 $4 
echo bfilterTile
nice -19 bin/bfilterTile $1 $2 $3 $4 $5 
echo bfilterTileOMP
nice -19 bin/bfilterTileOMP $1 $2 $3 $4 $5 
echo bfilterExpTile
nice -19 bin/bfilterExpTile $1 $2 $3 $4 $5 
echo bfilterExpTileOMP
bin/bfilterExpTileOMP $1 $2 $3 $4 $5 

