#!/usr/bin/env bash 

packinp=$(ls pack.inp 2> /dev/null | wc -l)
gaff=$(ls gaff.ff 2> /dev/null | wc -l)
polym=$(ls polym.in 2> /dev/null | wc -l)

[[ $packinp -eq 0 || $gaff -eq 0 || $polym -eq 0 ]] &&
echo "Make sure you have the following files in this dir:" &&
echo "- pack.inp" &&
echo "- gaff.ff" &&
echo "- polym.in" &&
echo "Make sure that pack.inp outputs a file named pack.xyz after" &&
echo "packmol has finished" &&
exit 1

# Need additional params for polymatic. The add_additional_params.py script is slow for a large system (uses itertools),
# so create a small system of 1 of each monomer, and add parameters to the small system instead.

mkdir small
cp pack.inp small/pack.inp
cp gaff.ff small/
# just 1 molecule for each different structure
$sed -i 's/number.*/number 1/' small/pack.inp
ls *xyz | xargs -I{} cp {} small/
cp polym.in small/
cd small
# clear output
{ packmol < pack.inp
lmp_gaff.py pack.xyz gaff.ff # creates pack.lmps
} > /dev/null
echo "Finding parameters that Polymatic needs. Could take about a minute..."
add_additional_params.py -l pack.lmps -f gaff.ff -p polym.in
# additional params now in small/pack.lmps
cd ..

## desired system ##
echo "Creating desired system"
{
packmol < pack.inp # generate pack.sh
lmp_gaff.py pack.xyz gaff.ff
} > /dev/null

add_correct_charges.py pack.lmps # Assigns charges using the GAMESS geodesic charge calculations

# Now take parameters from the small system and add into the desired datafile
echo "Adding additional parameters to datafile"
newtypes="$(grep 'types' small/pack.lmps)"
newparams="$(sed -n '/Masses/,/Atoms/p' small/pack.lmps | grep -v '^\s*Atoms')"
sed '/impropers/q' pack.lmps > tmp.start
grep 'lo.*hi' pack.lmps > tmp.boxsize
sed -n '/Atoms/,//p' pack.lmps > tmp.end

cat tmp.start <(echo "$newtypes") tmp.boxsize <(printf "\n$newparams\n\n") tmp.end > tmp.final
mv tmp.final pack.lmps

# Clean up
rm tmp* 
rm -rf small

polymatic_types.py pack.lmps > types.txt
# Molecule ID file needed for version of polymatic that considers free energy barriers between each monomer
generate_molecule_id_file.py
