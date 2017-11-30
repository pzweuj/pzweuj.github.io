for f in _data _includes _layouts _sass assets 404.html favicon.ico; do
    cp -Rfv ../jekyll-TeXt-theme/$f .;
done