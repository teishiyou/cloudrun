import os
from flask import Flask
from flask import send_from_directory

import io
from google.oauth2 import service_account
import googleapiclient.discovery
from googleapiclient.errors import HttpError
from googleapiclient.http import MediaIoBaseDownload
from apiclient.http import MediaFileUpload

app = Flask(__name__)


#sample to use it:
filename = 'test.csv'
filepath = '/app/benchling/screen/result/20230608/reads.xlsx'

folderid = '**'

@app.route("/")
def hello_world():
    SCOPES = ['https://www.googleapis.com/auth/drive']
    # サービスアカウントの鍵ファイル
    SERVICE_ACCOUNT_FILE = '/app/*.json'

    credentials = service_account.Credentials.from_service_account_file(
        SERVICE_ACCOUNT_FILE, scopes=SCOPES)

    delegated_credentials = credentials.with_subject('**')

    drive_service = googleapiclient.discovery.build('drive', 'v3', credentials=delegated_credentials)

    results = drive_service.files().list(
        corpora="drive",
        driveId='0ANctPunxD_ObUk9PVA',
        pageSize=10, 
        fields="files(id, name)",
        includeItemsFromAllDrives=True,
        supportsAllDrives=True
        ).execute()
    items = results.get('files', [])

    value = str(os.system("/usr/local/bin/Rscript /app/benchling/screen/run_screen.R"))
        
    file_metadata = {'name': filename,
                     "parents": [{'id': folderid}]
                     }
    media = MediaFileUpload(filepath, mimetype='text/xlsx',resumable=True)
    file = drive_service.files().create(body=file_metadata,
                                        media_body=media,
                                        fields='id').execute()

    return str(file.get('id'))

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))