import os

workingDirectory = '/home/rishubh/IITD/lab1'   #change this path
submissionDirectory = workingDirectory + '/submissionsl1p2'           # add your folder name instead of done

os.system('cat /dev/null > errorlog_part1.txt')
os.system('echo "entryno, test0, test1, test2, test3, test4, test5, test2_4, test3_4, test4_4, test5_4, test2_8, test4_8, test5_8, test2_16" > time_part1.csv')

errorlog = open('errorlog_part1.txt', 'w')

ctr = 0

for folder in os.listdir(submissionDirectory):
	ctr = ctr + 1

	studentDirectory = submissionDirectory + '/' + folder
	studentid =  folder

	errorlog.write('\n' + studentid)
	
	print (str(studentid) + ' started')

	contains = False
	for file in os.listdir(studentDirectory):
		if (file == 'convexhull.cpp'):
			contains = True

	if (contains == True):
		# os.system('cp ./A1P2_main.cpp ' + studentDirectory + '/')
		res_compile = os.system('g++ -fopenmp A1P2_main.cpp ' + studentDirectory + '/convexhull.cpp' + ' -o out1')
		
		if (res_compile == 0):
			exectimelog = str(studentid)

			errorflag = False

			res_exec = os.system('./out1 input0.txt 2 ' + studentDirectory + '/out2_0.txt  > ' + studentDirectory + '/res1')
			if (res_exec == 0):
				os.system('diff out0.txt ' + studentDirectory + '/out2_0.txt > diffres1')
				
				flag = False
				if (os.path.getsize('diffres1') > 0):
					flag = True
					errorflag = True
					errorlog.write(': answer mismatch (input0.txt : 2')

				resfile = open(studentDirectory + '/res1', 'r')
				t = resfile.read()
				ts = t.splitlines()
				if (flag) :
					exectimelog = exectimelog + (',' + str(25))
				else :
					exectimelog = exectimelog + (',' + str(ts[0]))
				resfile.close()
			else:
				errorflag = True
				# errorlog.write(': answer mismatch (' + testid + ' : ' + str(8) + ')')
				exectimelog = exectimelog + ',25'
			

            # run all test cases on 2 threads (8 test cases)
			for n in range(1,6):
				testid = 'input' + str(n) + '.txt'
				res_exec = os.system('./out1 ' + testid + ' 2 ' + studentDirectory + '/out2_' + str(n) + '.txt > ' + studentDirectory + '/res1')

				if (res_exec == 0):
					os.system('python l1p2.py out' + str(n) + '.txt ' + studentDirectory + '/out2_' + str(n) + '.txt > correct.txt')

					flag = False
					if (os.path.getsize('correct.txt') == 0):
						flag = True
						errorflag = True
						errorlog.write(': answer mismatch (' + testid + ' : ' + str(2) + ')')

					resfile = open(studentDirectory + '/res1', 'r')
					t = resfile.read()
					ts = t.splitlines()
					if (flag) :
						exectimelog = exectimelog + (',' + str(25))
					else :
						exectimelog = exectimelog + (',' + str(ts[0]))
					resfile.close()
				else:
					errorflag = True
					errorlog.write(': answer mismatch (' + testid + ' : ' + str(2) + ')')
					exectimelog = exectimelog + ',25'

			print('2 threads done')
			# on 4 thrds : 2,3,4,5
			for n in range(2, 6):
				testid = 'input' + str(n) + '.txt'
				res_exec = os.system('./out1 ' + testid + ' 4 ' + studentDirectory + '/out4_' + str(n) + '.txt > ' + studentDirectory + '/res1')

				if (res_exec == 0):
					os.system('python l1p2.py out' + str(n) + '.txt ' + studentDirectory + '/out4_' + str(n) + '.txt > correct.txt')

					flag = False
					if (os.path.getsize('correct.txt') == 0):
						flag = True
						errorflag = True
						errorlog.write(': answer mismatch (' + testid + ' : ' + str(4) + ')')

					resfile = open(studentDirectory + '/res1', 'r')
					t = resfile.read()
					ts = t.splitlines()
					if (flag) :
						exectimelog = exectimelog + (',' + str(25))
					else :
						exectimelog = exectimelog + (',' + str(ts[0]))
				else:
					errorflag = True
					errorlog.write(': answer mismatch (' + testid + ' : ' + str(4) + ')')
					exectimelog = exectimelog + ',25'

			# on 8 thrds : 2,4,5
			for n in [2,4,5]:
				testid = 'input' + str(n) + '.txt'
				res_exec = os.system('./out1 ' + testid + ' 8 ' + studentDirectory + '/out8_' + str(n) + '.txt > ' + studentDirectory + '/res1')

				if (res_exec == 0):
					os.system('python l1p2.py out' + str(n) + '.txt ' + studentDirectory + '/out8_' + str(n) + '.txt > correct.txt')

					flag = False
					if (os.path.getsize('correct.txt') == 0):
						flag = True
						errorflag = True
						errorlog.write(': answer mismatch (' + testid + ' : ' + str(8) + ')')

					resfile = open(studentDirectory + '/res1', 'r')
					t = resfile.read()
					ts = t.splitlines()
					if (flag) :
						exectimelog = exectimelog + (',' + str(25))
					else :
						exectimelog = exectimelog + (',' + str(ts[0]))
					resfile.close()
				else:
					errorflag = True
					errorlog.write(': answer mismatch (' + testid + ' : ' + str(8) + ')')
					exectimelog = exectimelog + ',25'

			#run test2 on 16thrds
			res_exec = os.system('./out1 input2.txt 16 ' + studentDirectory + '/out16_2.txt > ' + studentDirectory + '/res1')
			if (res_exec == 0):
				os.system('python l1p2.py out2.txt ' + studentDirectory + '/out16_2.txt > correct.txt')

				flag = False
				if (os.path.getsize('correct.txt') == 0):
					flag = True
					errorflag = True
					errorlog.write(': answer mismatch (input2.txt : 8)')

				resfile = open(studentDirectory + '/res1', 'r')
				t = resfile.read()
				ts = t.splitlines()
				if (flag) :
					exectimelog = exectimelog + (',' + str(25))
				else :
					exectimelog = exectimelog + (',' + str(ts[0]))
				resfile.close()
			else:
				errorflag = True
				errorlog.write(': answer mismatch (input2.txt : 8)')
				exectimelog = exectimelog + ',25'

			os.system('echo "' + exectimelog + '" >> time_part1.csv')

			if (errorflag == False):
				errorlog.write('pass')
		else:
			errorlog.write(": Compilation failed")
	else:
		errorlog.write(': Invalid name')
	
	print (str(studentid) + ' done')

errorlog.write("\n\nFolders read: " + str(ctr))
errorlog.close()

print ('Folders read:' + str(ctr))



