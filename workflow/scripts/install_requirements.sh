# sudo apt-get install unzip
conda_dir=$(conda info --base)

#git clone https://github.com/KennthShang/PhaBOX.git
#cd PhaBOX
#conda env create -f webserver.yml -n phabox

# database
fileid="1hjACPsIOqqcS5emGaduYvYrCzrIpt2_9"
filename="phagesuite_database.zip"
html=`curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid}"`
curl -Lb ./cookie "https://drive.google.com/uc?export=download&`echo ${html}|grep -Po '(confirm=[a-zA-Z0-9\-_]+)'`&id=${fileid}" -o ${filename}

# initial files
fileid="1E94ii3Q0O8ZBm7UsyDT_n06YekNtfV20"
filename="phagesuite_parameters.zip"
html=`curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid}"`
curl -Lb ./cookie "https://drive.google.com/uc?export=download&`echo ${html}|grep -Po '(confirm=[a-zA-Z0-9\-_]+)'`&id=${fileid}" -o ${filename}

unzip phagesuite_database.zip  > /dev/null
unzip phagesuite_parameters.zip  > /dev/null

# move the script to where the conda located
cp blastxml_to_tabular.py {conda_dir}/envs/phabox/bin/blastxml_to_tabular.py
chmod 777 {conda_dir}/envs/phabox/bin/blastxml_to_tabular.py