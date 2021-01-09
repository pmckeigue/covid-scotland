#!/bin/bash

pdftotext pharmaco.pdf
 
function wcpdfms() {
    local firstline=$( grep -n ^.\\?Background$ pharmaco.txt | sed 's/:.*//' )
    local lastline=$( grep -n ^.\\?Declarations$ pharmaco.txt | sed 's/:.*//' )
    local numlines=$(( lastline - firstline ))
    local wordcount=$( tail -n+$firstline pharmaco.txt | head -n $numlines | grep -v ^[0-9]\\+/[0-9]\\+$ | grep -v ^.20$ | wc -w )
    echo $wordcount
}

function wcabstr() {
    local firstline=$( grep -n "^Background .\\+" pharmaco.txt | sed 's/:.*//' )
    local lastline=$( grep -n "^.\\?Keywords" pharmaco.txt | sed 's/:.*//' )
    local numlines=$(( lastline - firstline ))
    local wordcount=$( tail -n+$firstline pharmaco.txt | head -n $numlines | grep -v ^[0-9]\\+/[0-9]\\+$ | grep -v ^.20$ | wc -w )
    echo $wordcount
}


echo "$(wcpdfms) words, abstract $(wcabstr) words"
