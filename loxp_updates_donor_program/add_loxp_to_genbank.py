from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
import os
def add_loxp(loxp_record):

    print("***********************shifting features for loxp**************************")

    #check for restriction sites and get location
    enzymes = [
        "BamHI",
        "HindIII",
        "EcoRI",
        "NotI",
    ]
    features = []
    
    restriction_enzyme_length = 0

    for feature in loxp_record.features:
        features.append(feature)
        if feature.type =="misc_feature" and str(feature.qualifiers.get('note')) in enzymes:
            restriction_site = feature
            restriction_enzyme_length = len(feature)

    lox_p_start_loc = restriction_site.location.end
    
    # Define loxP sequence and its location
    
    loxp_site_length = 34
    holomogy_arm_length = 40
    
    #gets end of the entire loxp feature including homology arms
    loxp_feature_end = lox_p_start_loc + loxp_site_length + holomogy_arm_length
    
    lox_p_strand = restriction_site.location.strand
    lox_p_seq =f"ATAACTTCGTATAGCATACATTATACGAAGTTAT"
    
    loxp_length = len(lox_p_seq) + restriction_enzyme_length
    
    lox_p_location_range = FeatureLocation(lox_p_start_loc, loxp_feature_end, strand=lox_p_strand)

    #move all the features down stream of loxp but keeps primers in the same spot
    for feature in features:
        if feature.location.start > loxp_feature_end and feature.type != "primer_bind":
            feature.location = FeatureLocation(feature.location.start + loxp_length, feature.location.end + loxp_length, strand=feature.location.strand)
    print(f"\n\n")
    
    #moves compound features like CDSs
    for feature in features:
        if isinstance(feature.location, CompoundLocation):
            updated_parts = []        
            for part in feature.location.parts:
                if loxp_feature_end > part.start and loxp_feature_end < part.end:
                    updated_parts.append(FeatureLocation(part.start, part.end + loxp_length, strand=part.strand))
                elif loxp_feature_end < part.start:
                    updated_parts.append(FeatureLocation(part.start + loxp_length, part.end + loxp_length, strand=part.strand))
                else:
                    updated_parts.append(part)
            updated_parts_sorted = sorted(updated_parts, key=lambda x: x.start)
    
            feature.location = CompoundLocation(updated_parts_sorted)
    
    print("loxP added successfully.\n\n")
    return loxp_record


if __name__ == "__main__":
    genbank = input("Enter the path to the GenBank file: ")
    guide_choice = input("Enter the guide name to add loxP after (e.g., 'g1'): ")
    add_loxp(genbank, guide_choice)
    
