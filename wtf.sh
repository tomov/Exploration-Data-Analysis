sbatch -p ncf --mem 500 -t 0-18:20 -o wtf_%j.out -e wtf_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'[m,i] = max([1,2,3]); exit'"
