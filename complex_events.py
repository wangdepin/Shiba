def complex_events(gtf_dic) -> list:
    """
    Make complex (CO) events list.
    These are complex events that don't fit into other categories
    and involve internal exons (not at transcript termini).

    Args:
        gtf_dic: A dictionary containing information about the GTF file.

    Returns:
        list: List of complex events
    """
    
    event_l = []
    for gene in gtf_dic.keys():
        if "intron_list" not in gtf_dic[gene]:
            continue
        chr = gtf_dic[gene]["chr"]
        strand = gtf_dic[gene]["strand"]
        gene_name = gtf_dic[gene]["gene_name"]
        intron_list = gtf_dic[gene]["intron_list"]
        intron_dic = gtf_dic[gene]["transcript_intron_dic"]
        exon_dic = gtf_dic[gene]["transcript_exon_dic"]
        
        # Compare all transcript pairs
        transcript_list = list(exon_dic.keys())
        for i in range(len(transcript_list)):
            for j in range(i+1, len(transcript_list)):
                tx1 = transcript_list[i]
                tx2 = transcript_list[j]
                
                # Get exons for each transcript
                tx1_exons = sorted(list(exon_dic[tx1]), key=lambda x: int(x.split(":")[1].split("-")[0]))
                tx2_exons = sorted(list(exon_dic[tx2]), key=lambda x: int(x.split(":")[1].split("-")[0]))
                
                # Skip if transcripts don't share any exons
                shared_exons = set(tx1_exons) & set(tx2_exons)
                if not shared_exons:
                    continue
                
                # Find differential regions between transcripts
                tx1_unique = list(set(tx1_exons) - set(tx2_exons))
                tx2_unique = list(set(tx2_exons) - set(tx1_exons))
                
                # Skip if no differences
                if not tx1_unique and not tx2_unique:
                    continue
                    
                # Skip if difference is only at transcript starts or ends
                # This will be handled by CF/CL functions
                
                # For CO, we need differences in internal exons
                # with shared flanking exons on both sides
                
                # Sort exons by position
                if tx1_unique:
                    tx1_unique = sorted(tx1_unique, key=lambda x: int(x.split(":")[1].split("-")[0]))
                if tx2_unique:
                    tx2_unique = sorted(tx2_unique, key=lambda x: int(x.split(":")[1].split("-")[0]))
                
                # Check if we have flanking exons
                # This is a simplified approach - in a full implementation
                # we would need more careful analysis to ensure these are truly internal
                if tx1_unique and tx2_unique:
                    tx1_start = min([int(x.split(":")[1].split("-")[0]) for x in tx1_unique])
                    tx1_end = max([int(x.split(":")[1].split("-")[1]) for x in tx1_unique])
                    tx2_start = min([int(x.split(":")[1].split("-")[0]) for x in tx2_unique])
                    tx2_end = max([int(x.split(":")[1].split("-")[1]) for x in tx2_unique])
                    
                    # Find shared exons before and after the unique regions
                    pre_common = []
                    post_common = []
                    
                    min_pos = min(tx1_start, tx2_start)
                    max_pos = max(tx1_end, tx2_end)
                    
                    for exon in shared_exons:
                        exon_end = int(exon.split(":")[1].split("-")[1])
                        exon_start = int(exon.split(":")[1].split("-")[0])
                        
                        if exon_end < min_pos:
                            pre_common.append(exon)
                        elif exon_start > max_pos:
                            post_common.append(exon)
                    
                    # If we have flanking exons on both sides, consider it a CO event
                    if pre_common and post_common:
                        # Determine which is included vs excluded based on exon count
                        if len(tx1_unique) > len(tx2_unique):
                            included_exons = tx1_unique
                            excluded_exons = tx2_unique
                            included_form_transcript = tx1
                            excluded_form_transcript = tx2
                        else:
                            included_exons = tx2_unique
                            excluded_exons = tx1_unique
                            included_form_transcript = tx2
                            excluded_form_transcript = tx1
                            
                        # Create event entry
                        exonlist_included = ";".join(included_exons)
                        exonlist_excluded = ";".join(excluded_exons)
                        pre_exon = max(pre_common, key=lambda x: int(x.split(":")[1].split("-")[1]))
                        post_exon = min(post_common, key=lambda x: int(x.split(":")[1].split("-")[0]))
                        
                        event_l.append([
                            exonlist_included, 
                            exonlist_excluded,
                            pre_exon,
                            post_exon,
                            strand, 
                            gene, 
                            gene_name,
                            included_form_transcript,
                            excluded_form_transcript
                        ])
    
    return event_l

def complex_first_events(gtf_dic) -> list:
    """
    Make complex first (CF) events list.
    These are complex events that involve the transcript start site.

    Args:
        gtf_dic: A dictionary containing information about the GTF file.

    Returns:
        list: List of complex first events
    """
    
    event_l = []
    for gene in gtf_dic.keys():
        if "intron_list" not in gtf_dic[gene]:
            continue
        chr = gtf_dic[gene]["chr"]
        strand = gtf_dic[gene]["strand"]
        gene_name = gtf_dic[gene]["gene_name"]
        intron_dic = gtf_dic[gene]["transcript_intron_dic"]
        exon_dic = gtf_dic[gene]["transcript_exon_dic"]
        
        # Compare all transcript pairs
        transcript_list = list(exon_dic.keys())
        for i in range(len(transcript_list)):
            for j in range(i+1, len(transcript_list)):
                tx1 = transcript_list[i]
                tx2 = transcript_list[j]
                
                # Get sorted exons for each transcript
                if strand == "+":
                    tx1_exons = sorted(list(exon_dic[tx1]), key=lambda x: int(x.split(":")[1].split("-")[0]))
                    tx2_exons = sorted(list(exon_dic[tx2]), key=lambda x: int(x.split(":")[1].split("-")[0]))
                    # Get first exons
                    tx1_first = tx1_exons[0] if tx1_exons else None
                    tx2_first = tx2_exons[0] if tx2_exons else None
                else:  # strand == "-"
                    tx1_exons = sorted(list(exon_dic[tx1]), key=lambda x: int(x.split(":")[1].split("-")[0]), reverse=True)
                    tx2_exons = sorted(list(exon_dic[tx2]), key=lambda x: int(x.split(":")[1].split("-")[0]), reverse=True)
                    # Get first exons
                    tx1_first = tx1_exons[0] if tx1_exons else None
                    tx2_first = tx2_exons[0] if tx2_exons else None
                
                # Skip if either transcript has no exons
                if not tx1_first or not tx2_first:
                    continue
                
                # Skip if first exons are the same
                if tx1_first == tx2_first:
                    continue
                
                # Find shared exons
                shared_exons = set(tx1_exons) & set(tx2_exons)
                if not shared_exons:
                    continue
                
                # Find the first shared exon - this will be the "post_common" exon
                if strand == "+":
                    shared_sorted = sorted(list(shared_exons), key=lambda x: int(x.split(":")[1].split("-")[0]))
                else:  # strand == "-"
                    shared_sorted = sorted(list(shared_exons), key=lambda x: int(x.split(":")[1].split("-")[0]), reverse=True)
                
                post_common = shared_sorted[0] if shared_sorted else None
                if not post_common:
                    continue
                
                # Get unique exons before the first shared exon
                tx1_unique = []
                tx2_unique = []
                
                for exon in tx1_exons:
                    if exon in shared_exons:
                        break
                    tx1_unique.append(exon)
                
                for exon in tx2_exons:
                    if exon in shared_exons:
                        break
                    tx2_unique.append(exon)
                
                # Skip if either transcript has no unique exons
                if not tx1_unique or not tx2_unique:
                    continue
                
                # Determine which is included vs excluded based on exon count
                if len(tx1_unique) > len(tx2_unique):
                    included_exons = tx1_unique
                    excluded_exons = tx2_unique
                    included_form_transcript = tx1
                    excluded_form_transcript = tx2
                else:
                    included_exons = tx2_unique
                    excluded_exons = tx1_unique
                    included_form_transcript = tx2
                    excluded_form_transcript = tx1
                    
                # Create event entry
                exonlist_included = ";".join(included_exons)
                exonlist_excluded = ";".join(excluded_exons)
                
                event_l.append([
                    exonlist_included, 
                    exonlist_excluded,
                    post_common,
                    strand, 
                    gene, 
                    gene_name,
                    included_form_transcript,
                    excluded_form_transcript
                ])
    
    return event_l

def complex_last_events(gtf_dic) -> list:
    """
    Make complex last (CL) events list.
    These are complex events that involve the transcript end site.

    Args:
        gtf_dic: A dictionary containing information about the GTF file.

    Returns:
        list: List of complex last events
    """
    
    event_l = []
    for gene in gtf_dic.keys():
        if "intron_list" not in gtf_dic[gene]:
            continue
        chr = gtf_dic[gene]["chr"]
        strand = gtf_dic[gene]["strand"]
        gene_name = gtf_dic[gene]["gene_name"]
        intron_dic = gtf_dic[gene]["transcript_intron_dic"]
        exon_dic = gtf_dic[gene]["transcript_exon_dic"]
        
        # Compare all transcript pairs
        transcript_list = list(exon_dic.keys())
        for i in range(len(transcript_list)):
            for j in range(i+1, len(transcript_list)):
                tx1 = transcript_list[i]
                tx2 = transcript_list[j]
                
                # Get sorted exons for each transcript
                if strand == "+":
                    tx1_exons = sorted(list(exon_dic[tx1]), key=lambda x: int(x.split(":")[1].split("-")[0]))
                    tx2_exons = sorted(list(exon_dic[tx2]), key=lambda x: int(x.split(":")[1].split("-")[0]))
                    # Get last exons
                    tx1_last = tx1_exons[-1] if tx1_exons else None
                    tx2_last = tx2_exons[-1] if tx2_exons else None
                else:  # strand == "-"
                    tx1_exons = sorted(list(exon_dic[tx1]), key=lambda x: int(x.split(":")[1].split("-")[0]), reverse=True)
                    tx2_exons = sorted(list(exon_dic[tx2]), key=lambda x: int(x.split(":")[1].split("-")[0]), reverse=True)
                    # Get last exons
                    tx1_last = tx1_exons[-1] if tx1_exons else None
                    tx2_last = tx2_exons[-1] if tx2_exons else None
                
                # Skip if either transcript has no exons
                if not tx1_last or not tx2_last:
                    continue
                
                # Skip if last exons are the same
                if tx1_last == tx2_last:
                    continue
                
                # Find shared exons
                shared_exons = set(tx1_exons) & set(tx2_exons)
                if not shared_exons:
                    continue
                
                # Find the last shared exon - this will be the "pre_common" exon
                if strand == "+":
                    shared_sorted = sorted(list(shared_exons), key=lambda x: int(x.split(":")[1].split("-")[0]), reverse=True)
                else:  # strand == "-"
                    shared_sorted = sorted(list(shared_exons), key=lambda x: int(x.split(":")[1].split("-")[0]))
                
                pre_common = shared_sorted[0] if shared_sorted else None
                if not pre_common:
                    continue
                
                # Get unique exons after the last shared exon
                tx1_unique = []
                tx2_unique = []
                
                found_shared = False
                for exon in reversed(tx1_exons):
                    if not found_shared and exon in shared_exons:
                        found_shared = True
                        continue
                    if found_shared:
                        tx1_unique.append(exon)
                
                found_shared = False
                for exon in reversed(tx2_exons):
                    if not found_shared and exon in shared_exons:
                        found_shared = True
                        continue
                    if found_shared:
                        tx2_unique.append(exon)
                
                # Skip if either transcript has no unique exons
                if not tx1_unique or not tx2_unique:
                    continue
                
                # Determine which is included vs excluded based on exon count
                if len(tx1_unique) > len(tx2_unique):
                    included_exons = tx1_unique
                    excluded_exons = tx2_unique
                    included_form_transcript = tx1
                    excluded_form_transcript = tx2
                else:
                    included_exons = tx2_unique
                    excluded_exons = tx1_unique
                    included_form_transcript = tx2
                    excluded_form_transcript = tx1
                    
                # Create event entry
                exonlist_included = ";".join(included_exons)
                exonlist_excluded = ";".join(excluded_exons)
                
                event_l.append([
                    exonlist_included, 
                    exonlist_excluded,
                    pre_common,
                    strand, 
                    gene, 
                    gene_name,
                    included_form_transcript,
                    excluded_form_transcript
                ])
    
    return event_l