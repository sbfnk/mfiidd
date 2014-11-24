#!/bin/sh

CURRENT_DIR=$(pwd)

for I in $*;
do

    if [ ! -f $I ]; then
        echo "File $I not found"
        continue
    fi

    filename=$(basename $I)
    filebase=$(basename $I .md)
    filepath=$(dirname $I)

    cd $filepath

    # create pdf version
    pandoc --smart -f markdown-markdown_in_html_blocks   -o $filebase.pdf $filename
    # create HTML version
    pandoc --toc  --mathjax  -c github-pandoc.css -f markdown -t html -o $filebase.html $filename

    # convert links to HTML
    sed -i "" 's/\.md/.html/g' $filebase.html

    cd $CURRENT_DIR

done
