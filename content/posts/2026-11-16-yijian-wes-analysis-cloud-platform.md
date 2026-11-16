---
title: "YiJian: A Cloud WES Analysis Platform for Germline Research"
tags: software
---

YiJian is a cloud WES analysis platform for germline research, combining a FASTQ-to-results pipeline, quality control, variant annotation, and browser-based result review.

Whole-exome sequencing generates data. Turning those data into a usable set of research candidates requires an analysis workflow, reference resources, annotation tools, and a way to inspect the results. For a small team, operating that infrastructure can become a project of its own.

I built **YiJian** to bring these steps into one web interface. Researchers can use the hosted platform without deploying the analysis stack themselves. Teams that need to operate their own infrastructure can explore the open-source, self-hosted edition.

- **[Explore the YiJian cloud platform](https://yijian.schema-bio.com)**
- **[Read the English user guide](https://www.schema-bio.com/en/docs/yijian/user-guide)**
- **[Read the self-hosted deployment guide](https://www.schema-bio.com/en/docs/yijian/deployment)**

The hosted service is for **Research Use Only (RUO)**. It provides bioinformatics analysis and research results, rather than clinical diagnostic reports or individual medical advice. Service availability and applicable terms should be checked before submitting data.

## A WES analysis pipeline, with a platform around it

A **WES analysis pipeline** runs the computational steps needed to process whole-exome sequencing data. A **WES analysis platform** also needs to help people submit data, follow tasks, inspect quality metrics, review variants, and retrieve results.

YiJian brings those activities together. Its workflow starts with sample registration and sequencing data, continues through analysis and annotation, and ends with results that researchers can browse and export.

The practical goal is to reduce the manual handoffs between scripts, spreadsheets, job logs, and separate viewers. Users still need to evaluate data quality and biological evidence, but they can do that within a connected workspace.

## From FASTQ to annotated germline variants

### Prepare samples and sequencing data

Users register samples, associate their sequencing files, and create analysis tasks through the web interface. The data center provides a place to manage inputs, while the task center shows analysis progress and results.

In a self-hosted installation, teams can also register files already available on mounted storage. That option is useful for large datasets that do not need to be uploaded again through a browser.

![Sample registration in the YiJian WES analysis platform](https://www.schema-bio.com/images/docs/yijian/sample-info.png)

### Run the whole-exome sequencing analysis workflow

YiJian uses WDL-based bioinformatics workflows, with containerized tools and reference resources configured for the selected analysis. The underlying stack includes tools such as GATK, DeepVariant, and VEP; the tools and outputs used for a particular task depend on its workflow configuration.

During development, I moved the execution architecture toward CPU-based workflows. The reason was operational: I wanted an analysis system that I could deploy and maintain without depending on the availability of a particular GPU instance.

I also moved the workflow definitions to WDL and used miniwdl for execution. Separating the web interface, scheduler, and compute agent made it possible to manage tasks without exposing a management port on every compute node.

These choices underpin the platform. They are not a promise of a fixed turnaround time: runtime depends on sequencing volume, workflow settings, compute resources, and storage performance.

### Check quality and review variant annotations

After a task completes, users can inspect quality-control results and browse annotated variants. Gene lists and filters help narrow a broad result set to candidates relevant to the research question.

The broader YiJian software stack includes result views for SNVs and small insertions/deletions, CNVs, STRs, mobile element insertions, mitochondrial variants, UPD, and runs of homozygosity. Available result categories depend on the selected workflow, its inputs, and the service configuration.

Support for a result category does not imply equivalent detection performance across all loci or datasets. Capture design, coverage, reference resources, and method limitations remain part of the evaluation.

![SNV and indel annotation results in YiJian](https://www.schema-bio.com/images/docs/yijian/snp-indel.png)

### Inspect reads with an embedded IGV viewer

An annotation table is a starting point for review. Alignment context can help researchers investigate coverage, mismatches, and other signals around a candidate variant.

YiJian integrates an IGV viewer so users can inspect reads alongside the result interface, when the relevant alignment and index files are available. This connects candidate filtering with visual examination of the underlying sequencing evidence.

![Embedded IGV viewer for inspecting sequencing reads around a candidate variant](https://www.schema-bio.com/images/docs/yijian/igv.png)

## Cloud WES analysis or self-hosted deployment?

The hosted and self-hosted editions offer different operating arrangements. Choosing between them starts with how your team wants to manage infrastructure and data.

| Consideration | Hosted cloud platform | Open-source self-hosted edition |
| --- | --- | --- |
| Infrastructure | Use the hosted analysis environment | Deploy and maintain your own environment |
| Getting started | Check service eligibility, create an account, and prepare authorized research data | Prepare compute, storage, references, and service configuration |
| Data retention | Follow the hosted service's retention and deletion policies | Establish retention and backup policies for your deployment |
| Intended output | Research analysis results; no clinical diagnostic reports | Configure workflows and integrations for your independently governed use case |

The open-source edition includes family-analysis workflows and configurable report-service integration. Those capabilities should not be read as a promise that every feature is offered by the hosted service. Self-hosting also does not, by itself, establish regulatory compliance or make an output a validated diagnostic result.

For teams evaluating **cloud-based exome analysis**, the hosted platform offers a way to use the workflow without first installing the full stack. For teams that need control over their environment, the deployment guide explains the alternative.

## Data handling is part of the workflow

The hosted platform is an analysis service, not a permanent sequencing-data archive. Its current policy defines retention periods for samples, tasks, input files, and results. Users should review the [YiJian privacy policy](https://yijian.schema-bio.com/privacy) and service terms, and download the materials they need within the applicable windows.

Do not include direct personal identifiers in sample names, filenames, or free-text fields. Using a sample code does not necessarily make genetic data anonymous. Teams must have the authority to submit their data and assess the processing arrangements that apply to their research.

AI-assisted features, where enabled, help organize evidence and identify research candidates. Their output requires independent review and must not be treated as a diagnosis or a substitute for professional judgment.

## Open-source components behind YiJian

At the end of 2025, I planned to build both a germline exome platform and a tumor analysis platform in my spare time. As the work progressed, I narrowed the scope to germline exome analysis. Completing one usable system became more valuable than maintaining two unfinished ones.

The resulting stack is modular, and its core components are publicly available:

| Component | Purpose | License |
| --- | --- | --- |
| [YiJian](https://github.com/SchemaBio/YiJian) | Web interface for analysis tasks and result review | Apache-2.0 |
| [Octopus](https://github.com/SchemaBio/Octopus) | Backend scheduling and task management | MIT |
| [schema-germline](https://github.com/SchemaBio/schema-germline) | WDL-based germline analysis workflows | MIT |
| [Sepiida](https://github.com/SchemaBio/Sepiida) | Compute-node status collection and reporting | MIT |
| [SepiLens](https://github.com/SchemaBio/SepiLens) | AI-assisted evidence retrieval and assessment | Apache-2.0 |

Publishing the components gives technical teams a way to inspect the implementation, understand the workflow, and contribute improvements. The self-hosted edition also provides a path for users whose infrastructure requirements differ from those of the hosted service.

## Getting started with YiJian

If you are evaluating a WES analysis platform, start with one concrete question: can it take your input data through the workflow you need and produce outputs you can independently assess?

Read the [YiJian user guide](https://www.schema-bio.com/en/docs/yijian/user-guide) to understand the sample, task, and result workflow. Then [visit the cloud platform](https://yijian.schema-bio.com) to review access and service information, or follow the [self-hosted deployment guide](https://www.schema-bio.com/en/docs/yijian/deployment) to run the open-source stack on your own infrastructure.

YiJian connects the computational pipeline with the everyday work around it: preparing data, tracking analysis, checking quality, and reviewing germline variants. That is the platform I set out to build, and the foundation I will continue to improve.
