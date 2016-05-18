cd structures/
find . -name '*' | xargs rm 
cd ../unafold_files/
find . -name '*' | xargs rm 
cd ..
rm -rf transterm_files/*
