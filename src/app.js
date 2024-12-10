"use strict";
let page = 0;
const style = document.createElement("style");
style.textContent = site.f["style.css"];
document.head.appendChild(style);
var update_time = document.timeline.currentTime;
document.body.insertAdjacentHTML("beforeend", site.f["page.html"]);
const rest = document.getElementById("rest");
const hash = decodeURI(globalThis.location.hash.substring(1)) || "";
let new_page = hash == "sectionzero" ? 0 : page;
rest.insertAdjacentHTML("afterbegin", '<div id="page0"></div>\n');
for (let i = page - 1; i >= 0; i--) {
  rest.insertAdjacentHTML("afterbegin", `<div id="page${i}"></div>\n`);
  if (hash != "" && site.slugs[`page${i}`].content.includes(`id="${hash}"`)) {
    new_page = i;
  }
}
for (let i = page + 1; i <= 39; i++) {
  rest.insertAdjacentHTML("beforeend", `<div id="page${parseInt(i)}"></div>\n`);
  if (hash != "" && site.slugs[`page${i}`].content.includes(`id="${hash}"`)) {
    new_page = i;
  }
}
page = new_page;
const scroll_nav = (entries) => {
  const section = entries[0].target.id;
  document.getElementById("TOC").querySelectorAll("a").forEach((a) => {
    if (a.id.slice(4) == section) {
      if (a.parentElement.parentElement.parentElement.localName == "details") {
        a.parentElement.parentElement.parentElement.setAttribute("open", "");
      }
      a.style.fontSize = "2em";
    } else {
      a.style.fontSize = ".9em";
    }
    a.scrollIntoView({
      block: "center",
    });
  });
};
const options = { rootMargin: "-50% 0% -50% 0%" };
const observer = new IntersectionObserver(scroll_nav, options);
const s_dest = page == 0 && hash == ""
  ? "sectionzero"
  : hash == ""
  ? `page${page}`
  : hash;
refresh_page();
document.getElementById(s_dest).scrollIntoView();
const pd = document.getElementById("pandoc");
pd.addEventListener("scroll", () => {
  if (document.timeline.currentTime - update_time > 200) {
    const el_loc = document.getElementById("page" + parseInt(page))
      .getBoundingClientRect();
    if (el_loc.top > pd.clientHeight || el_loc.bottom < 0) {
      if (el_loc.bottom < 0) {
        page = page + 1;
        refresh_page();
      } else {
        page = page - 1;
        refresh_page();
      }
    }
  }
});
function refresh_page() {
  for (let i = -2; i <= 2; i++) {
    const pe = document.getElementById(`page${page + i}`);
    if (site.slugs[`page${page + i}`] && pe.innerHTML.length == 0) {
      pe.innerHTML = `${site.slugs[`page${page + i}`].content}\n</div>\n`;
      pe.querySelectorAll("h1,h2,h3").forEach(
        (heading_link) => {
          observer.observe(heading_link);
        },
      );
    }
  }
  if (page <= 2) {
    page_0_top();
    top_page(true);
  } else {
    top_page(false);
  }
}
function top_page(state) {
  if (state) {
    for (const i of ["logo", "titlehead", "subtitle"]) {
      document.getElementById(i).style.display = "flex";
    }
    for (const i of ["abstract", "thanks", "barbottom"]) {
      document.getElementById(i).style.display = "block";
    }
    document.getElementById("pdf").style.display = "flex";
    document.getElementById("term_legal").addEventListener("click", () => {
      document.getElementById("full").style.display = "none";
      document.getElementById("legal").style.display = "block";

      document.getElementById("return").addEventListener("click", () => {
        document.getElementById("full").style.display = "grid";
        document.getElementById("legal").style.display = "none";
      });
    });
  } else {
    for (const i of ["logo", "titlehead", "subtitle"]) {
      document.getElementById(i).style.display = "none";
    }
    for (const i of ["pdf", "abstract", "thanks", "barbottom"]) {
      document.getElementById(i).style.display = "none";
    }
  }
}
globalThis.addEventListener("hashchange", () => {
  const hash = decodeURI(globalThis.location.hash.substring(1)) || "sectionzero";
  let found_page = 0;
  for (const i in site.slugs) {
    const page_int = parseInt(i.slice(4));
    if (site.slugs[i].content.includes('id="' + hash + '"')) {
      found_page = page_int;
    }
  }
  page = found_page;
  refresh_page();
  document.getElementById(hash).scrollIntoView();
  return false;
});
document.getElementById("sean_button").addEventListener("click", () => {
  page = 0;
  refresh_page();
  history.pushState({}, globalThis.title, "#sectionzero");
  document.getElementById("sectionzero").scrollIntoView();
});
document.getElementById("mini_nav").innerHTML = site.f["mininav.html"];
document.getElementById("TOC").innerHTML = site.f["toc.html"];
document.getElementById("legal").innerHTML = site.f["legal.html"];
function page_0_top() {
  if (document.getElementById("abstract").innerHTML.length == 0) {
    document.getElementById("abstract").innerHTML = site.f["abstract.html"];
    document.getElementById("thanks").innerHTML = site.f["thanks.html"];
    document.getElementById("top_graph").innerHTML = site.f["topgraph.html"];
    document.getElementById("logo").innerHTML = site.f["logo.html"];
    document.getElementById("titlehead").innerHTML = site.f["titlehead.html"];
    document.getElementById("pdf").innerHTML = site.f["pdf.html"];
    document.getElementById("subtitle").innerHTML = site.f["subtitle.html"];
  }
}
