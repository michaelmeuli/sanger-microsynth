"""Fetch Microsynth result mails over IMAP and save their attachments."""
from __future__ import annotations

import email
import imaplib
import json
import re
import zipfile
from email.message import Message
from email.utils import parsedate_to_datetime
from pathlib import Path

from . import MAIL_DIR, STATE_FILE

SEQ_SUFFIXES = {".ab1", ".seq", ".fasta", ".fa", ".fna", ".txt", ".pdf", ".zip"}


def _load_state() -> dict:
    if STATE_FILE.exists():
        return json.loads(STATE_FILE.read_text())
    return {}


def _save_state(state: dict) -> None:
    STATE_FILE.parent.mkdir(parents=True, exist_ok=True)
    STATE_FILE.write_text(json.dumps(state, indent=2))


def _safe_name(name: str) -> str:
    return re.sub(r"[^\w.\-]+", "_", name).strip("._") or "attachment"


def _safe_extract(zip_path: Path, dest: Path) -> list[Path]:
    """Extract only files whose resolved path stays inside dest (zip-slip guard)."""
    out = []
    with zipfile.ZipFile(zip_path) as zf:
        for info in zf.infolist():
            if info.is_dir():
                continue
            target = (dest / info.filename).resolve()
            if dest.resolve() not in target.parents:
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(zf.read(info))
            out.append(target)
    return out


def _save_attachments(msg: Message, dest: Path) -> list[Path]:
    saved = []
    for part in msg.walk():
        filename = part.get_filename()
        if not filename or part.get_content_maintype() == "multipart":
            continue
        name = _safe_name(filename)
        if Path(name).suffix.lower() not in SEQ_SUFFIXES:
            continue
        dest.mkdir(parents=True, exist_ok=True)
        path = dest / name
        path.write_bytes(part.get_payload(decode=True) or b"")
        saved.append(path)
        if path.suffix.lower() == ".zip":
            try:
                saved.extend(_safe_extract(path, dest))
            except zipfile.BadZipFile:
                pass
    return saved


def fetch_new(cfg: dict[str, str]) -> list[Path]:
    """Download attachments of not-yet-processed mails whose From header
    contains SENDER_FILTER. Uses BODY.PEEK so the mailbox is left unchanged.
    Returns the new batch directories (one per mail with attachments)."""
    state = _load_state()
    new_dirs: list[Path] = []
    with imaplib.IMAP4_SSL(cfg["IMAP_HOST"]) as imap:
        imap.login(cfg["IMAP_USER"], cfg["IMAP_PASSWORD"])
        imap.select(cfg["IMAP_FOLDER"], readonly=True)
        _, data = imap.search(None, "FROM", f'"{cfg["SENDER_FILTER"]}"')
        for num in data[0].split():
            _, fetched = imap.fetch(num, "(BODY.PEEK[])")
            msg = email.message_from_bytes(fetched[0][1])
            msg_id = (msg.get("Message-ID") or "").strip()
            if not msg_id or msg_id in state:
                continue
            try:
                date = parsedate_to_datetime(msg["Date"]).strftime("%Y%m%d")
            except (TypeError, ValueError):
                date = "undated"
            batch = MAIL_DIR / f"{date}_{_safe_name(msg_id)[:40]}"
            files = _save_attachments(msg, batch)
            state[msg_id] = {"subject": str(msg.get("Subject", "")), "date": date,
                             "dir": batch.name, "files": [f.name for f in files]}
            if any(f.suffix.lower() == ".ab1" for f in files):
                new_dirs.append(batch)
    _save_state(state)
    return new_dirs
