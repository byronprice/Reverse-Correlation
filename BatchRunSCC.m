% BatchRunSCC.m

myCluster = parcluster('local');

if getenv('ENVIRONMENT')
   myCluster.JobStorageLocation = getenv('TMPDIR'); 
end

parpool(myCluster,8);

Animals = [34271,34272,43650,43652,43653,62500,62501,62502];

parfor ii=1:8
    RevCorrMov(Animals(ii),20170421,'pink',1);
end

delete(gcp);