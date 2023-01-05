function writeModel(A,b)

    dlmwrite('A.txt',A,'\t');
    dlmwrite('b.txt',b,'\n');

end