import mysql.connector


##### User data #####

# MySQL Connection Information
Database = 'PFAMphylostratigraphy'
User = 'anhnguyenphung'
Host = '127.0.0.1'
Password = '9400_changethis'

# table to get sequences from
AlignmentTable = 'PfamAlignments'
PfamColumn = "PfamUID"
UIDColumn = "UID" ### pfams can be present multiple times in a genome

Seq = 'AlignedPeptide'


# The script then establishes a connection with the SQL database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)

mycursor = cnx.cursor(buffered = True)


pfam = 'PF00450'
  
    
PullProteinsStatement = "SELECT * FROM %s WHERE %s='%s'" %(AlignmentTable, PfamColumn, pfam)
mycursor.execute(PullProteinsStatement)
ProteinSequences = mycursor.fetchall()

print(ProteinSequences)


