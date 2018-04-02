import requests
import time ### TIME STOP
import sys
import textwrap


def pairwise_align_needle_api(s1, s2, smatrix, dna):

    mtx = define_matrix(smatrix)

    GAP_COST = 100
    GAP_EXT = 0.0005
    END_WEGT = "false"
    END_OPEN = 1
    TYPE = "dna" if dna  else "protein"

    url = 'https://www.ebi.ac.uk/Tools/services/rest/emboss_needle/run/'
    url_status = 'https://www.ebi.ac.uk/Tools/services/rest/emboss_needle/status/'
    url_result = 'https://www.ebi.ac.uk/Tools/services/rest/emboss_needle/result/'

    headers = {'Host': 'www.ebi.ac.uk', 'Content-Type':'application/x-www-form-urlencoded', 'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/60.0.3112.113 Safari/537.36'}
    data_info = {'bsequence' : s1,
                 'email' : 'bio01@gmail.com',
                 'stype': TYPE,
                 'asequence':s2,
                 'matrix': mtx,
                 'gapopen' : GAP_COST
                }

    try:


        sys.stdout.write("Conectando API....\n\n")

        req = requests.post(url, headers=headers, params=data_info)

        req.raise_for_status() ## Erro se o status code fro diferente de 200

        job_id = req.text
        url_status = url_status + job_id

        try:
            while True:
                req = requests.get(url_status, headers=headers)
                req.raise_for_status() ## Erro se o status code fro diferente de 200

                if req.text == 'FINISHED':
                    break

                time.sleep(2)

            url_result = url_result + job_id + '/aln'

            try:
                req = requests.get(url_result, headers=headers)
                req.raise_for_status() ## Erro se o status code fro diferente de 200
                resultado_final = req.text


                resultado_final_lines = resultado_final.splitlines()

                resultado = {}
                i=1
                resultado['Seq_1'] = ''
                resultado['Seq_2'] = ''

                for line in resultado_final_lines:
                    if "# Length:" in line:
                        resultado['Lenght'] = int(line[9:])
                    if "# Identity:" in line:
                        resultado['Identity'] = line[11:].strip()
                    if "# Similarity:  " in line:
                        resultado['Similarity'] = line[15:].strip()
                    if "# Gaps:  " in line:
                        resultado['Gaps'] = line[9:].strip()
                    if "# Score: " in line:
                        resultado['Score'] = int(float(line[8:]))
                    if line.startswith("EMBOSS_001"):
                        if (i % 2 == 0 ):
                            resultado['Seq_1'] =resultado['Seq_1'] + (line[21:]).split("     ")[0]
                        else:
                            resultado['Seq_2'] =resultado['Seq_2'] + (line[21:]).split("     ")[0]
                        i += 1


                print("--- EMBOSS API---")
                print("Gap cost {}".format(GAP_COST))
                print("Score Matrix {}".format(smatrix))
                print("Score {}".format(resultado['Score']))
                #print(resultado['Seq_1'])
                #print(resultado['Seq_2'])
                result1 = ''.join(i for i in resultado['Seq_1'] if not i.isdigit())
                result1 = result1.strip()
                result1 = result1.replace(" ", "")

                result2 = ''.join(i for i in resultado['Seq_2'] if not i.isdigit())
                result2 = result2.strip()
                result2 = result2.replace(" ", "")


                textwrap.wrap(result1, width=100)
                print(textwrap.fill(result1, width=100))
                textwrap.wrap(result2, width=100)
                print(textwrap.fill(result2, width=100))

            except requests.exceptions.HTTPError as err:
                print(err)
        except requests.exceptions.HTTPError as err:
            print(err)
    except requests.exceptions.HTTPError as err:
        print(err)


def define_matrix(m):
    if (m == 'pam30'):
        return "EPAM30"
    elif (m == 'blosum62'):
        return "EBLOSUM62"
    elif (m == 'pam250'):
        return "EPAM250"
    elif (m == 'dnafull'):
        return "EDNAFULL"
    else:
        return "EBLOSUM50"








