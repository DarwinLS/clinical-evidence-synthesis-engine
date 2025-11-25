from fastapi import FastAPI, Request, Form
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
import uvicorn

app = FastAPI()

# Telling FastAPI where the CSS and HTML files are
app.mount("/static", StaticFiles(directory="static"), name="static")
templates = Jinja2Templates(directory="templates")

# The Home Page (GET Request)
@app.get("/", response_class=HTMLResponse)
async def read_root(request: Request):
    return templates.TemplateResponse("index.html", {"request": request, "result": None})

# The Search Action (POST Request)
@app.post("/analyze", response_class=HTMLResponse)
async def analyze_supplement(
    request: Request, 
    supplement: str = Form(...), 
    age: int = Form(...)
):
    # MOCK PIPELINE (Simulates what the "Writer LLM" will eventually return)
    mock_result = {
        "supplement": supplement,
        "demographic": f"Age {age}",
        "summary": [
            {
                "topic": "Cognitive Performance",
                "text": "Evidence suggests a significant benefit for short-term memory tasks.",
                "citation_id": "12345",
                "citation_meta": "RCT, n=45, 2023"
            },
            {
                "topic": "Safety Profile",
                "text": "No adverse effects noted in the standard dosage range.",
                "citation_id": "67890",
                "citation_meta": "Meta-Analysis, n=1200, 2021"
            }
        ]
    }

    # Send mock data back to the HTML page
    return templates.TemplateResponse("index.html", {
        "request": request, 
        "result": mock_result,
        "search_term": supplement
    })

# Entry point for running the app directly
if __name__ == "__main__":
    uvicorn.run("main:app", host="127.0.0.1", port=8000, reload=True)