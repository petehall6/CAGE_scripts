a = "blarh blah blah. More blah.  Foo.  Bar. Lead-Pete Hall."
c = "Lead-Pete Hall."


print(f"a:  {len(a.split('.'))}")

lst = [a,c]

for let in lst:
    lead = let.split("Lead-")[1]
    print(lead)