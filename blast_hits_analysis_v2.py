
from datetime import datetime, date

result1_string = open("blast hits results/f1.txt").read().splitlines()
result2_string = open("blast hits results/f2.txt").read().splitlines()
result3_string = open("blast hits results/f3.txt").read().splitlines()
result4_string = open("blast hits results/f4.txt").read().splitlines()
result5_string = open("blast hits results/f5.txt").read().splitlines()
result6_string = open("blast hits results/f6.txt").read().splitlines()

result1_string = ["frame1 \t" + x for x in result1_string]
result2_string = ["frame2 \t" + x for x in result2_string]
result3_string = ["frame3 \t" + x for x in result3_string]
result4_string = ["frame4 \t" + x for x in result4_string]
result5_string = ["frame5 \t" + x for x in result5_string]
result6_string = ["frame6 \t" + x for x in result6_string]

all_result = result1_string + result2_string  + result3_string  + result4_string  + result5_string  + result6_string

best_result = []
#best_result.append("")
temp = [""]
i = 1
for z in range(0,len(all_result)):
    item = all_result[z] 
    temp = item.split("\t")

    if i==1:
        frame = temp[0]
        title = temp[1]
        accession = temp[2]
        evalue = temp[11]
        bitsocre = temp[12]
        best_result.append( "frame: " + frame +"title: " +title + " -accession: " + accession + " -E Value: " + evalue + " -Bit score " + bitsocre )

    else:
  
        if title == temp[1]:
            if float(evalue) > float(temp[11]):
                frame = temp[0]
                title = temp[1]
                accession = temp[2]
                evalue = temp[11]
                bitsocre = temp[12]
                best_result.pop()
                best_result.append("frame: " + frame + "title: " +title + " -accession: " + accession + " -E Value: " + evalue + " -Bit score " + bitsocre )
        else:
            frame = temp[0]
            title = temp[1]
            accession = temp[2]
            evalue = temp[11]
            bitsocre = temp[12]
            status= "notfound"
            for items in best_result:
                line = items.split(" ")
                past_title = line[3]
                if title == past_title:
                    status ="found"
                    break
            
            if status=="notfound":
                best_result.append( "frame: " + frame + "title: " +title + " -accession: " + accession + " -E Value: " + evalue + " -Bit score " + bitsocre )

            

    for y in range(z,len(all_result)):
        temp_item = all_result[y] 
        temp_line = temp_item.split("\t")
        temp_title = temp_line[1]
        temp_evalue = temp_line[11]
        if title in temp_line:
            if float(evalue) > float(temp_evalue):
                best_result.pop()
                frame = temp_line[0]
                title = temp_line[1]
                accession = temp_line[2]
                evalue = temp_line[11]
                bitsocre = temp_line[12]
                best_result.append( "frame: " + frame + "title: " +title + " -accession: " + accession + " -E Value: " + evalue + " -Bit score " + bitsocre )



    
    i = i + 1


best_result.insert(0,"")

current_dateTime = datetime.now()

x=current_dateTime
f = open("blast analysis/mahmoodvand_6frame_translation_result" + str(x.year) +"_"+ str(x.month)+"_"+str(x.day)+"_"+str(x.hour)+"_"+str(x.minute) + ".txt", "w")
f.write("the result of best blast hit from contig by Mohamadreza Mahmoodvand - " + str(current_dateTime) )

for line in best_result:
    f.write(line)
    f.write('\n')
f.close()

print(best_result)

print("end")