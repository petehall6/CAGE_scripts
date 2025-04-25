import subprocess

#*  BE SURE TO TYPE conda activate cage into the terminal before running!!! *#


subprocess.run('pip install pandas')

subprocess.run('pip install numpy')

subprocess.run('pip install -r requirements.txt')

subprocess.run('pip install pywin32')

subprocess.run('pip install ttkbootstrap')


print("Package installation complete.")