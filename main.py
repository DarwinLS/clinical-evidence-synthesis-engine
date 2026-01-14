from fastapi import FastAPI, Request, Form
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
import uvicorn
import logging

# IMPORT PIPELINE MODULES
from src.services.pubmed import search_pubmed
from src.pipeline.extractor import extract_metadata
from src.pipeline.scorer import score_studies
from src.pipeline.synthesizer import synthesize_report

# Configure logging so you can see what's happening in the terminal
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI()

app.mount("/static", StaticFiles(directory="static"), name="static")
templates = Jinja2Templates(directory="templates")

@app.get("/", response_class=HTMLResponse)
async def read_root(request: Request):
    return templates.TemplateResponse("index.html", {"request": request, "result": None})

@app.post("/analyze", response_class=HTMLResponse)
async def analyze_supplement(
    request: Request, 
    supplement: str = Form(...), 
    age: int = Form(...),
    goal: str = Form("general") # Default to 'general' if not provided
):
    logger.info(f"--- STARTING PIPELINE for {supplement} (Age {age}, Goal: {goal}) ---")
    
    # 1. FETCH
    # We fetch more results (e.g., 20) to ensure we have enough after filtering
    raw_studies = search_pubmed(supplement, max_results=20)
    
    if not raw_studies:
        return render_error(request, "No clinical studies found on PubMed for this query.")

    logger.info(f"Step 1 Complete: Fetched {len(raw_studies)} raw studies.")

    # 2. EXTRACT
    # Worker LLM parses text -> JSON
    extracted_data = extract_metadata(raw_studies)
    
    if not extracted_data:
        return render_error(request, "Failed to extract metadata from studies.")

    logger.info(f"Step 2 Complete: Extracted metadata.")

    # 3. SCORE & FILTER
    # Logic Engine picks winners based on Age + Quality
    top_studies = score_studies(extracted_data, user_age=age, goal=goal)
    
    if not top_studies:
        return render_error(request, "Studies found, but none matched your age group/quality criteria.")

    logger.info(f"Step 3 Complete: Selected Top {len(top_studies)} studies.")

    # 4. SYNTHESIZE
    final_report = synthesize_report(supplement, age, top_studies, goal=goal)
    
    if not final_report:
        return render_error(request, "Failed to generate summary report.")

    logger.info("Step 4 Complete: Report generated successfully.")

    # SANITIZE CITATIONS
    # Create the lookup dictionary first
    study_lookup = {s['id']: s for s in top_studies}

    # Clean LLM output before sending to Frontend
    if final_report and "summary" in final_report:
        for section in final_report["summary"]:
            clean_ids = []
            raw_ids = section.get("citation_ids", [])
            
            # Handle case where LLM returns a single string instead of a list
            if isinstance(raw_ids, str):
                raw_ids = [raw_ids]
            
            for raw_id in raw_ids:
                # Force string to match dictionary keys
                str_id = str(raw_id)
                
                # Only keep the ID if we actually have the study metadata
                if str_id in study_lookup:
                    clean_ids.append(str_id)
            
            # Update the section with the clean list
            section["citation_ids"] = clean_ids

    return templates.TemplateResponse("index.html", {
        "request": request, 
        "result": final_report,
        "study_lookup": study_lookup,
        "search_term": supplement,
        "search_age": age,
        "search_goal": goal
    })

def render_error(request, message):
    """Helper to show error messages on the frontend"""
    return templates.TemplateResponse("index.html", {
        "request": request, 
        "error": message
    })

if __name__ == "__main__":
    uvicorn.run("main:app", host="127.0.0.1", port=8000, reload=True)