# get the minimal representation for one allele
minrep_onerow = function(pos, ref, alt) {
    ref = unlist(strsplit(ref,split="")) # string to vector
    alt = unlist(strsplit(alt,split="")) # string to vector
    if (length(ref) == 1 && length(alt) == 1) {
        # do nothing
    } else {
        # strip off identical suffixes
        while(alt[length(alt)] == ref[length(ref)] && min(length(alt),length(ref)) > 1) {
            alt = alt[1:length(alt)-1]
            ref = ref[1:length(ref)-1]
        }
        # strip off identical prefixes and increment position
        while(alt[1] == ref[1] && min(length(alt),length(ref)) > 1) {
            alt = alt[2:length(alt)]
            ref = ref[2:length(ref)]
            pos = pos + 1
        }
    }
    ref = paste(ref,collapse="")
    alt = paste(alt,collapse="")
    return (list(pos, ref, alt))
}
 
# vectorized version of above
minrep_vectorized = function(pos, ref, alt) {
    stopifnot(length(pos) == length(ref) && length(ref) == length(alt))
    for (i in 1:length(pos)) {
        list_of_results = minrep_onerow(pos[i],ref[i],alt[i])
        pos[i] = list_of_results[[1]]
        ref[i] = list_of_results[[2]]
        alt[i] = list_of_results[[3]]
        # if (i %% 1000 == 0) {
        #   print(i)
        # }
    }
    return (data.frame(pos, ref, alt))
}