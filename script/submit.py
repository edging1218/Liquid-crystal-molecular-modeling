import subprocess
import os
dic = {'n': 78, 
       'l': 77, 
       'r':25}
dir_name = 'n' + str(dic['n']) + '_l' + str(dic['l']) + '_r' + str(dic['r'])
cur_dir = os.path.dirname(os.path.realpath(__file__))
subprocess.call('mkdir '+dir_name, shell=True)
os.chdir(os.path.join(cur_dir, dir_name))
with open('../param_template') as f:
    template = f.read()
with open('param.in', 'w') as f:
    f.write(template % dic)

with open('../sbatch_template') as f:
    template = f.read()
with open('depablo.sbatch', 'w') as f:
    f.write(template % dic)
# subprocess.call('sbatch depablo.sbatch', shell=True)
os.chdir(cur_dir)
