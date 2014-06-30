#!/bin/sh

if [ ! -f $1 ]; then
    echo "File $1 not found"
    exit 1
fi

CURRENT_DIR=$(pwd)

filename=$(basename $1)
filebase=$(basename $1 .md)
filepath=$(dirname $1)

cd $filepath

# create pdf version
pandoc --smart -f markdown-markdown_in_html_blocks   -o $filebase.pdf $filename
# create HTML version
pandoc --toc  --mathjax  -c github-pandoc.css -f markdown -t html -o $filebase.html $filename

# convert links to HTML
sed -i "" 's/\.md/.html/g' $filebase.html

cd $CURRENT_DIR
