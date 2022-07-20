function numw = get_parallel_numworkers()

if verLessThan('matlab', '8.3')
    numw =  matlabpool('size') ;
else
    p=gcp('nocreate');
    if isempty(p);
        numw=0;
    else
        numw = p.NumWorkers;
    end
end
