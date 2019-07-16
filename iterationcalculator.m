function iterationcalculator(coloriter)
if(ischar(coloriter))
    %decompose col
    colvec=double(coloriter);
    colvec=colvec-96;
    msk = colvec<0;
    colvec(msk)=0;
    colnum = colvec(1).*(26.^2) + colvec(2).*(26.^1) + colvec(3).*(26.^0);
    
    fprintf("Iteration Number: %4.0f\n",colnum);
    eachdiv = 10000/16;
    workernumber=ceil(colnum/eachdiv)+1;
    fprintf("Worker Number: %2.0f\n",workernumber);
else
    coloriter=coloriter+1;
    out(1) = floor(coloriter./(26.^2));
    sub = out(1)*(26^2);
    newnum=coloriter-sub;
    out(2) = floor(newnum./(26));
    sub = out(2)*26;
    out(3) = newnum-sub;
    msk = out==0;
    out(msk)=32-96;
    out=char(out+96);
    fprintf("Column: %s\n",out);
    
end

end
