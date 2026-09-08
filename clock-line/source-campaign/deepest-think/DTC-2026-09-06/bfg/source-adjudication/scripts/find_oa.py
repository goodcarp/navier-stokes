import json,urllib.request,urllib.parse,sys
dois=["10.1088/1361-6544/ac62de","10.1088/1361-6544/ad68b9","10.1007/s00021-025-00947-x","10.1007/s40818-025-00199-y","10.1007/s00021-024-00888-x"]
for d in dois:
    u="https://api.openalex.org/works/doi:"+d+"?mailto=contact@example.invalid"
    j=json.load(urllib.request.urlopen(u))
    print(j['publication_year'],'|',j['title'])
    for l in j.get('locations',[]):
        if l.get('is_oa'): print('   OA:',l.get('pdf_url') or l.get('landing_page_url'))
