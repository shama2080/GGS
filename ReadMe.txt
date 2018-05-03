I want to use ggs algorithm in a C# project, So i implement it in c# visual studio 2017.
To write your own C# function that uses ggs,
	1. Add Accord library to your references.
	2. create a new object of cGGS (MyGGSClass) and then add the following code to your function:
		MyGGSClass.fGGS(out breaks,out plotPoints, data, m_kMax, m_dLamb, nFeats);



* Some unused functions may be incomplete like fMulti_run_wrapper.