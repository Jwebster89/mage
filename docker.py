import subprocess
import sys, os
from gooey import Gooey
from gooey import GooeyParser
import pathlib


class Mage():
	def __init__(self,path,csv):
		self.path = path
		self.csv = csv

	def batch_creation(self, files, csv):
		(shortname, ext) =os.path.splitext(csv)
		run_id=os.path.basename(shortname)
		with open(os.path.join(pathlib.Path(__file__).parent.absolute(),"DOUBLE_CLICK_ME_mage.bat"), "w") as bat_file:
			bat_file.write("docker volume create mage_data\n")
			bat_file.write("xcopy /E/H/C/I {} \\\wsl$\docker-desktop-data\\version-pack-data\community\docker\\volumes\mage_data\_data \n".format(files))
			bat_file.write("docker run -it -v mage_data:/home/data/files --name jwebster89/mage mage --path /home/data/files --csv /home/data/files/{} --force\n".format(os.path.basename(csv)))
			bat_file.write("mkdir ./mage_output_{} \n".format(run_id))
			bat_file.write("xcopy /E/H/C/I \\\wsl$\docker-desktop-data\\version-pack-data\community\docker\\volumes\mage_data\_data .\mage_output_{} \n".format(run_id))
			bat_file.write("docker rm -f mage\n")
			bat_file.write("docker volume rm -f mage_data\n")

	def bat_run(self):
		bat_file=os.path.join(pathlib.Path(__file__).parent.absolute(),"mage.bat")
		subprocess.run([bat_file])

	def run(self):
		self.batch_creation(self.path, self.csv)
		#self.bat_run()

@Gooey
def main():
	parser = GooeyParser(description="MLSA Generator")
	parser.add_argument('--path', help="path to ab1 and ref files", widget='DirChooser')
	parser.add_argument('--csv', help="csv file with sequence information", widget='FileChooser') 

	args=parser.parse_args()
	path=os.path.normpath(args.path)
	csv_input=args.csv

	job=Mage(path, csv_input)
	job.run()

if __name__ == '__main__':
	main()
