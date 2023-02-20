from Bio.SeqUtils import six_frame_translations
from datetime import datetime, date


import Bio
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

def blast(protein_seq,evalue):
    result_handle = NCBIWWW.qblast("blastp", "nr", protein_seq)
    from Bio.Blast import NCBIXML
    blast_records = NCBIXML.parse(result_handle)
    blast_records_list = list(blast_records)
    E_VALUE_THRESH =float(evalue)
    count = 0
    count_obj = 0
    current_dateTime = datetime.now()
    x=current_dateTime

    #f = open("mahmoodvand_result" + str(x.year) +"_"+ str(x.month)+"_"+str(x.day)+"_"+str(x.hour)+"_"+str(x.minute) + ".txt", "w")
    #f.write("the result of contig blast by Mohamadreza Mahmoodvand - " + str(current_dateTime) )
    final_result =[]

    for obj in blast_records_list:
        print()
        final_result.append("")
        print("-------------------------start of protein examination-------------------------------")
        final_result.append("-------------------------start of protein examination-------------------------------")
        print("protein seq: "+ str(count_obj+1), obj.query)
        final_result.append("protein seq: "+ str(count_obj+1) + str(obj.query))
        count = 0
        for alignment in obj.alignments:
            if count<6:
                for hsp in alignment.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        print()
                        final_result.append("")
                        print("****Protein Matched: " + str(count+1) +"****")
                        final_result.append("****Protein Matched: " + str(count+1) +"****")
                        print("title:", alignment.title)
                        final_result.append("title:" + str(alignment.title))
                        print("bit score:", hsp.bits)
                        final_result.append("bit score:" + str( hsp.bits))
                        print("total score:", hsp.score)
                        final_result.append("total score:" + str( hsp.score))
                        print("E value:", hsp.expect)
                        final_result.append("E value:"+ str(hsp.expect))
                        print()
                        final_result.append("")
                        count = count +1
            else:
                break
        if(count==0):
            print(" Not Found any Protein")
            final_result.append(" Not Found any Protein")
        count_obj = count_obj + 1
        print("-------------------------end of protein examination-------------------------------")
        final_result.append("-------------------------end of protein examination-------------------------------")
        print()
        final_result.append("")

    return final_result

def pro_maker(six_frame):
    result = six_frame.splitlines()

    del result[0:4]

    line_counts = len(result)

    spilted_seq = line_counts // 10

    pro_seq = ["","","","","",""]

    p1_index = 1
    p2_index = 2
    p3_index = 3
    p4_index = 6
    p5_index = 7
    p6_index = 8

    for i in range(0,spilted_seq):
        pro_seq[0] =  pro_seq[0].strip() + result[p1_index].strip() 
        pro_seq[1] = pro_seq[1].strip() + result[p2_index].strip() 
        pro_seq[2] =  pro_seq[2].strip() + result[p3_index].strip() 
        pro_seq[3] =  pro_seq[3].strip() + result[p4_index].strip() 
        pro_seq[4] = pro_seq[4].strip() + result[p5_index].strip() 
        pro_seq[5] =  pro_seq[5].strip() + result[p6_index].strip()

        p1_index = p1_index + 10
        p2_index = p2_index + 10
        p3_index = p3_index + 10
        p4_index = p4_index + 10
        p5_index = p5_index + 10
        p6_index = p6_index + 10

   

    return pro_seq

# Main Program

try:
    current_dateTime = datetime.now()

    x=current_dateTime

    f = open("6frame/mahmoodvand_6frame_translation_result" + str(x.year) +"_"+ str(x.month)+"_"+str(x.day)+"_"+str(x.hour)+"_"+str(x.minute) + ".txt", "w")
    f.write("the result of 6 Frame extraction from contig by Mohamadreza Mahmoodvand - " + str(current_dateTime) )
    final_result =[]

    fasta_string = open("100contig").read().splitlines()

    pro_candid = []

    z=1
    c_count = 1
    blast_result = []
    print("Hi, This is Contig Parser by Mohamadreza Mahmoodvand")
    print("How many contig do you want to parse?")
    stop_limit = int(input("type your stop limit (number) :  "))
    E_Value = float(input("type your considered E Value :  "))
    blast_option = int(input("Do you want Blast 6 frames (type 1) or not (type 0) : "))
    format_file = int(input("select format file: for fastA type 1 : "))

    final_result.append("")

    if(format_file == 1):
        final_result1 =[]
        final_result2 =[]
        final_result3 =[]
        final_result4 =[]
        final_result5 =[]
        final_result6 =[]

        for z in range(1,len(fasta_string),2):
            dna_seq =fasta_string[z]
                       
            pro_candid = pro_maker(six_frame_translations(dna_seq))
            Pro_count =1
            for item in pro_candid:
                if(Pro_count == 1):
                    final_result1.append(fasta_string[z-1])
                    final_result1.append(item)
                if(Pro_count == 2):
                    final_result2.append(fasta_string[z-1])
                    final_result2.append(item)
                if(Pro_count == 3):
                    final_result3.append(fasta_string[z-1])
                    final_result3.append(item)
                if(Pro_count == 4):
                    final_result4.append(fasta_string[z-1])
                    final_result4.append(item)
                if(Pro_count == 5):
                    final_result5.append(fasta_string[z-1])
                    final_result5.append(item)
                if(Pro_count == 6):
                    final_result6.append(fasta_string[z-1])
                    final_result6.append(item)

                Pro_count = Pro_count + 1

        f1 = open("fastA/mahmoodvand_fastA_1_" + str(x.year) +"_"+ str(x.month)+"_"+str(x.day)+"_"+str(x.hour)+"_"+str(x.minute) + ".fa", "w")
        f2 = open("fastA/mahmoodvand_fastA_2_" + str(x.year) +"_"+ str(x.month)+"_"+str(x.day)+"_"+str(x.hour)+"_"+str(x.minute) + ".fa", "w")
        f3 = open("fastA/mahmoodvand_fastA_3_" + str(x.year) +"_"+ str(x.month)+"_"+str(x.day)+"_"+str(x.hour)+"_"+str(x.minute) + ".fa", "w")
        f4 = open("fastA/mahmoodvand_fastA_4_" + str(x.year) +"_"+ str(x.month)+"_"+str(x.day)+"_"+str(x.hour)+"_"+str(x.minute) + ".fa", "w")
        f5 = open("fastA/mahmoodvand_fastA_5_" + str(x.year) +"_"+ str(x.month)+"_"+str(x.day)+"_"+str(x.hour)+"_"+str(x.minute) + ".fa", "w")
        f6 = open("fastA/mahmoodvand_fastA_6_" + str(x.year) +"_"+ str(x.month)+"_"+str(x.day)+"_"+str(x.hour)+"_"+str(x.minute) + ".fa", "w")
    

        for line1 in final_result1:
            f1.write(line1)
            f1.write('\n')      
        for line2 in final_result2:
            f2.write(line2)
            f2.write('\n')
        for line3 in final_result3:
            f3.write(line3)
            f3.write('\n')  
        for line4 in final_result4:
            f4.write(line4)
            f4.write('\n')  
        for line5 in final_result5:
            f5.write(line5)
            f5.write('\n')  
        for line6 in final_result6:
            f6.write(line6)
            f6.write('\n')        

        f1.close()
        f2.close()
        f3.close()
        f4.close()
        f5.close()
        f6.close()  

    else:
        for z in range(1,len(fasta_string),2):
            dna_seq =fasta_string[z]
            print("<<<------------------START of Blast of Protein seq------------------>>>")
            final_result.append("")
            final_result.append("<<<-----------------START of Blast of Protein seq------------------>>>")
            print( str(c_count) + "- Contig Name :" + fasta_string[z-1]  )
            final_result.append(str(c_count) + "- Contig Name :" + fasta_string[z-1])
            pro_candid = pro_maker(six_frame_translations(dna_seq))
            print("candidate 6 frame proteins are ...")
            final_result.append("candidate 6 frame proteins are ...")
            p_count = 1
            for item in pro_candid:
                print( str(p_count) + " : " + item)
                final_result.append(str(p_count) + " : " + item)
                p_count = p_count + 1

                if(blast_option == 1):
                    try:
                        blast_result = blast(item,E_Value)
                    except:
                        blast_result = ["Network Error - Connection Failed due to Iran internet Access - NO BLAST for protein !"]
                    for item2 in blast_result:
                        final_result.append(item2)
                        print(item2)
                    

                
            print("<<<-----------------END of Blast of Protein seq------------------->>>")
            final_result.append("<<<-----------------END of Blast of Protein seq------------------->>>")
            
            c_count = c_count + 1

            if((z // 2) > stop_limit):
                break
        

    #result  = pro_maker(six_frame_translations("GAGATAGTATTTATGGAATATGATTGTATTTTTATAGAACGTTATCACCTATATGTATCCTTGTTCTTCAGCGTCAGCTCCATGTCAACTTAATCAAGGTTCATCTGCAATTCCTAGTGCCTAGTCTCATCGGCTGTTAATCGGTCCTACTCCAGCTTAGCTTTCTTCACATATATTTACTCTATACGCATATGATCCCAGAATCCCCATGGCACTTATATAATCTGGTTCATTCTCTTCAACAACGATTTCTACAGTTGCTTTTACGACAGTTTTTGAAACAGTATTTGCTATCTCTTTTTCATCATCGTCATAACCATAAATAAAAACGTTCTCATCCTTTGATGATTCCATCCCACCGCTTGGCCCGGTAACAGTAAAGTTCTCCACATCACCGTTTTCATCAAGTTCTATTAAGTAGTCACATGCTTGATGCGCCACATCCATCTCTAAAGCCGAACTACCCAAACTGATTCCTACCGCTTTATCATAT"))


    #print(result)

    end_dateTime = datetime.now()

    print("Start time :",str(current_dateTime))
    print("End time :",str(end_dateTime))
    final_result.append("Number of Translated contig: " + str(c_count -1))
    final_result.append("Start time :"+str(current_dateTime))
    final_result.append("End time :"+str(end_dateTime))
    final_result.append("E Value: " + str(0.05))
    dur =  end_dateTime - current_dateTime
    final_result.append("Duration :"+str(dur))
    print("duration: "+ str(dur))
    final_result.append("End of file")

    for line in final_result:
        f.write(line)
        f.write('\n')

    # for item2 in blast_result:
    #     f.write(item2)
    #     print(item2)
    #     f.write('\n')


    #f.writelines(final_result)
    f.close()
except:
    for line in final_result:
        f.write(line)
        f.write('\n')
print("enjoy! you can find the report at 6frame Folder")
print("end of the program")