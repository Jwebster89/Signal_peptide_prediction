#!/usr/bin/env python3

import os, sys, csv
import argparse
import subprocess
import matplotlib.pyplot as plt
from venn import venn

class Signal_Peptides():
	def __init__(self,protein,phobius,signalp,deepsig,wolfpsort):
		self.protein=protein
		self.phobius=phobius
		self.signalp=signalp
		self.deepsig=deepsig
		self.wolfpsort=wolfpsort

	def result_dir(self):
		if not os.path.exists("Signal_Peptide_Prediction"):
			os.mkdir("Signal_Peptide_Prediction")

	def run_phobius(self,protein):
		(protein_prefix,ext)=os.path.splitext(protein)
		base=os.path.basename(protein_prefix)
		if os.path.exists(f"Signal_Peptide_Prediction/Phobius/{base}.phobius") and self.phobius:
			print(f"Output for {protein} and Phobius already exists. Using previous output.")
			return(os.path.join("Signal_Peptide_Prediction/Phobius/",base+".phobius"))
		else:
			if self.phobius:
				if not os.path.exists("Signal_Peptide_Prediction/Phobius/"):
					os.mkdir("Signal_Peptide_Prediction/Phobius")
				print(f"Currently Running Phobius on {protein}")
				with open(os.path.join("Signal_Peptide_Prediction/Phobius/"+base+".phobius"), 'w') as out_phobius_h:
					subprocess.run(['./bin/tmp/tmpe49ws9/phobius/phobius.pl','-short',protein], stdout=out_phobius_h, stderr=subprocess.DEVNULL)
				return(os.path.join("Signal_Peptide_Prediction/Phobius/",base+".phobius"))

	def run_signalp(self, protein):
		(protein_prefix,ext)=os.path.splitext(protein)
		base=os.path.basename(protein_prefix)
		if os.path.exists(f"Signal_Peptide_Prediction/SignalP/{base}_summary.signalp5") and self.signalp:
			print(f"Output for {protein} and SignalP already exists. Using previous output.")
			return(str(os.path.join("Signal_Peptide_Prediction/SignalP/",base+"_summary.signalp5")))
		else:
			if self.signalp:
				if not os.path.exists("Signal_Peptide_Prediction/SignalP"):
					os.mkdir("Signal_Peptide_Prediction/SignalP")
				print(f"Currently Running SignalP on {protein}")
				subprocess.run(['signalp','-fasta',protein,'-org','euk','-prefix',"Signal_Peptide_Prediction/SignalP/"+base], stdout=subprocess.DEVNULL)
				return(str(os.path.join("Signal_Peptide_Prediction/SignalP/",base+"_summary.signalp5")))

	def run_deepsig(self, protein):
		(protein_prefix,ext)=os.path.splitext(protein)
		base=os.path.basename(protein_prefix)
		if os.path.exists(f"Signal_Peptide_Prediction/Deepsig/{base}.deepsig") and self.deepsig:
			print(f"Output for {protein} and Deepsig already exists. Using previous output")
			return(os.path.join("Signal_Peptide_Prediction/Deepsig/",base+".deepsig"))
		else:
			if self.deepsig:
				if not os.path.exists("Signal_Peptide_Prediction/Deepsig"):
					os.mkdir("Signal_Peptide_Prediction/Deepsig")
				print(f"Currently Running Deepsig on {protein}")
				(protein_prefix,ext)=os.path.splitext(protein)
				base=os.path.basename(protein_prefix)
				subprocess.run(['deepsig','-f',protein, '-o',"Signal_Peptide_Prediction/Deepsig/"+base+".deepsig", '-k','euk'], stdout=subprocess.DEVNULL)
				return(os.path.join("Signal_Peptide_Prediction/Deepsig/",base+".deepsig"))

	def run_wolfpsort(self, protein):
		(protein_prefix,ext)=os.path.splitext(protein)
		base=os.path.basename(protein_prefix)
		if os.path.exists(f"Signal_Peptide_Prediction/WolfPSort/{base}.wolfpsort") and self.wolfpsort:
			print(f"Output for {protein} WolfPSort already exists. Using previous output")
			return(os.path.join("Signal_Peptide_Prediction/WolfPSort/",base+".wolfpsort"))
		else:
			if self.wolfpsort:
				if not os.path.exists("Signal_Peptide_Prediction/WolfPSort"):
					os.mkdir("Signal_Peptide_Prediction/WolfPSort")
				print(f"Currently Running WolfPSort on {protein}")
				(protein_prefix,ext)=os.path.splitext(protein)
				base=os.path.basename(protein_prefix)
				input=subprocess.Popen(['cat',protein],stdout=subprocess.PIPE)
				with open(os.path.join("Signal_Peptide_Prediction/WolfPSort/"+base+".wolfpsort"), 'w') as out_wolfpsort_h:
					subprocess.run(['./bin/WoLFPSort/bin/runWolfPsortSummary','fungi'],stdin=input.stdout,stdout=out_wolfpsort_h)
				return(os.path.join("Signal_Peptide_Prediction/WolfPSort/",base+".wolfpsort"))


	def read_phobius(self, input):
		if self.phobius:
			print("reading Phobius")
			with open(input, 'r') as input_p:
				phobius_list=[]
				rd = csv.reader(input_p, delimiter=" ", quotechar='"')
				for row in rd:
					if row[-2] == "Y" and row[-4] == "0":
						phobius_list.append(row[0])
			return(phobius_list, "phobius")

	def read_signalp(self,input):
		if self.signalp:
			print("reading signalP")
			with open(input, 'r') as input_sp:
				signalp_list=[]
				rd=csv.reader(input_sp, delimiter="\t", quotechar='"')
				for row in rd:
					if "SP" in row[1]:
						signalp_list.append(row[0])
			return(signalp_list, "signalp")

	def read_deepsig(self, input):
		if self.deepsig:
			print("reading deepsig")
			with open(input, 'r') as input_ds:
				deepsig_list=[]
				rd=csv.reader(input_ds, delimiter="\t", quotechar='"')
				for row in rd:
					if row[2] == "Signal peptide":
						deepsig_list.append(row[0])
			return(deepsig_list, "deepsig")
	
	def read_wolfpsort(self, input):
		if self.wolfpsort:
			print("reading WolfPSort")
			with open(input, 'r') as input_wps:
				wps_list=[]
				rd=csv.reader(input_wps, delimiter=" ", quotechar='"')
				for row in rd:
					if "extr" in row[1]:
						wps_list.append(row[0])
			return(wps_list,"wolfpsort")
	
	def venn(self, protein, *args):
		if not os.path.exists("Signal_Peptide_Prediction/Venn"):
			os.mkdir("Signal_Peptide_Prediction/Venn")
		set_list=[]
		set_names=[]
		for program in args:
			if program is None:
				continue
			else:
				(list,name)=program
				set_list.append(set(list))
				set_names.append(name)
		colours=["tab:blue","tab:orange","tab:green","tab:purple"]
		dict={}
		for i in range(len(set_names)):
			dict[set_names[i]] = set_list[i]
		venn(dict,cmap=colours)
		(protein_prefix,ext)=os.path.splitext(protein)
		base=os.path.basename(protein_prefix)
		plt.title(f'{base} Signal Peptide Prediction')
		plt.savefig(f'Signal_Peptide_Prediction/Venn/{base}.venn.png')
		print("Done!")
		return(set_list,base)
	
	def output_intersection(self,prefix, sets):
		inter={}
		for j, item in enumerate(sets):
			if j == 0:
				inter=item
			elif j== len(sets):
				break
			else:
				inter=inter.intersection(item)
		with open(f"Signal_Peptide_Prediction/{prefix}.prediction_consensus.txt", 'w') as f_out:
			f_out.write('\n'.join(inter))
		print("Done!")

	def run(self):
		print(f"\n=== Running Prediction on {self.protein} ===")
		self.result_dir()
		phobius_out=self.run_phobius(self.protein)
		signalp_out=self.run_signalp(self.protein)
		deepsig_out=self.run_deepsig(self.protein)
		wolfpsort_out=self.run_wolfpsort(self.protein)
		print("=== Reading output of Predictions ===")
		listp=self.read_phobius(phobius_out)
		listsp=self.read_signalp(signalp_out)
		listds=self.read_deepsig(deepsig_out)
		listwps=self.read_wolfpsort(wolfpsort_out)
		count=0
		for val in [listp, listsp,listds,listwps]:
			if val !=None:
				count+=1
		if count > 1:
			print("=== Writing Venn Diagram ===")
			(sets,prefix)=self.venn(self.protein, listp,listsp, listds, listwps)
			print("=== Writing consensus list ===")
			self.output_intersection(prefix, sets)


def main():
	parser = argparse.ArgumentParser(description="Signal Peptide Prediction", add_help=False)
	required = parser.add_argument_group('Required Arguments')
	required.add_argument('-i','--input',type=str,required=True)
	
	optional = parser.add_argument_group('Optional Arguments')
	optional.add_argument('-p', '--phobius', required=False, action='store_true', help="Run Phobius")
	optional.add_argument('-s', '--signalp', required=False, action='store_true', help="Run SignalP")
	optional.add_argument('-d', '--deepsig', required=False, action='store_true', help="Run Deepsig")
	optional.add_argument('-w', '--wolfpsort', required=False, action='store_true', help="Run WolfPSort")
	optional.add_argument("-h", "--help", action="help", help="show this help message and exit")

	args=parser.parse_args()
	input=args.input
	phobius=args.phobius
	signalp=args.signalp
	deepsig=args.deepsig
	wolfpsort=args.wolfpsort
	job = Signal_Peptides(input,phobius,signalp,deepsig,wolfpsort)
	job.run()

if __name__ == '__main__':
	main()
