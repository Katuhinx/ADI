connect:
	ssh mic@dvmh.keldysh.ru
connect2:
	ssh student_3@alex-freenas.ddns.net -p2000
send:
	scp -r -P 2000 adi_utils *.c *.h *.cdv Makefile student_3@alex-freenas.ddns.net:/home/student_3/adi
send_mic:
	scp *.c *.h *.cdv mic@dvmh.keldysh.ru:/home/dvmuser4/Katrine/adi
clear:
	rm *.o
compile_init:
	gcc -o adi adi.c adi_utils/adi_utils.c
compile_oacc:
	ompcc -o himenoOACC himenoOACC.c
compile_xmp:
	xmpcc -o himenoXMP himenoXMP.c
compile_xacc:
	xmpcc -xacc -o himenoXACC himenoXACC.c

compile_dvmh:
	./dvm c -o adiD adiD.cdv adi_utils/adi_utils.c

load:
	scp -P 2000 student_3@alex-freenas.ddns.net:/home/student_3/himeno/*.c /mnt/d/himeno
	scp -P 2000 student_3@alex-freenas.ddns.net:/home/student_3/himeno/*.h /mnt/d/himeno


