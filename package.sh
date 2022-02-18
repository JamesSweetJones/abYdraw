version=V1.0
tmpdir=abYdraw_${version}
mkdir $tmpdir
cp abYdraw.py install.py $tmpdir
tar zcvf ${tmpdir}.tgz ${tmpdir}
\rm -rf ${tmpdir}
