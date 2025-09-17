# n8n Workflow Templates for GenomeAMRAnalyzer
# Ready-to-import workflow definitions

## 1. Core Analysis Storage Workflow

```json
{
  "name": "AMR Analysis to Database",
  "nodes": [
    {
      "parameters": {
        "path": "amr-analysis-complete",
        "options": {}
      },
      "name": "Analysis Complete Webhook",
      "type": "n8n-nodes-base.webhook",
      "typeVersion": 1,
      "position": [240, 300],
      "webhookId": "amr-analysis-complete"
    },
    {
      "parameters": {
        "functionCode": "// Parse and validate AMR analysis results\nconst analysis = items[0].json;\n\n// Validate required fields\nif (!analysis.metadata || !analysis.summary) {\n  throw new Error('Invalid analysis data structure');\n}\n\n// Extract key information\nconst processedData = {\n  analysis_id: analysis.metadata.analysis_id,\n  timestamp: analysis.metadata.timestamp,\n  user_email: analysis.metadata.user_email,\n  genome_count: analysis.summary.total_genomes,\n  resistance_genes: analysis.summary.total_resistance_genes,\n  high_priority_findings: analysis.summary.high_priority_mutations,\n  quality_score: analysis.quality_metrics?.overall_score || 0.0,\n  processing_time: analysis.summary.processing_time_minutes,\n  files: analysis.files,\n  findings: analysis.findings\n};\n\n// Check for high-priority alerts\nprocessedData.needs_alert = processedData.high_priority_findings >= 5;\n\nreturn [{json: processedData}];"
      },
      "name": "Parse Analysis Data",
      "type": "n8n-nodes-base.function",
      "typeVersion": 1,
      "position": [460, 300]
    },
    {
      "parameters": {
        "operation": "insert",
        "table": "amr_analyses",
        "columns": "analysis_id,user_email,timestamp,genome_count,resistance_genes,high_priority_findings,quality_score,processing_time",
        "values": "={{$json.analysis_id}},={{$json.user_email}},={{$json.timestamp}},={{$json.genome_count}},={{$json.resistance_genes}},={{$json.high_priority_findings}},={{$json.quality_score}},={{$json.processing_time}}"
      },
      "name": "Store to Database",
      "type": "n8n-nodes-base.postgres",
      "typeVersion": 1,
      "position": [680, 300],
      "credentials": {
        "postgres": {
          "id": "amr_database",
          "name": "AMR Results Database"
        }
      }
    },
    {
      "parameters": {
        "conditions": {
          "boolean": [
            {
              "value1": "={{$json.needs_alert}}",
              "value2": true
            }
          ]
        }
      },
      "name": "Check Alert Needed",
      "type": "n8n-nodes-base.if",
      "typeVersion": 1,
      "position": [900, 300]
    },
    {
      "parameters": {
        "subject": "🚨 High-Priority AMR Findings - Analysis {{$json.analysis_id}}",
        "message": "Analysis {{$json.analysis_id}} has detected {{$json.high_priority_findings}} high-priority antimicrobial resistance findings.\n\nSummary:\n- Genomes Analyzed: {{$json.genome_count}}\n- Resistance Genes Found: {{$json.resistance_genes}}\n- Quality Score: {{$json.quality_score}}\n- Processing Time: {{$json.processing_time}} minutes\n\nReview the detailed results in your analysis dashboard.",
        "options": {}
      },
      "name": "Send Alert Email",
      "type": "n8n-nodes-base.emailSend",
      "typeVersion": 2,
      "position": [1120, 200],
      "credentials": {
        "smtp": {
          "id": "smtp_config",
          "name": "SMTP Configuration"
        }
      }
    },
    {
      "parameters": {
        "operation": "create",
        "bucket": "amr-analysis-results",
        "fileName": "={{$json.analysis_id}}/analysis_summary.json",
        "fileContent": "={{JSON.stringify($json, null, 2)}}",
        "options": {}
      },
      "name": "Archive to S3",
      "type": "n8n-nodes-base.awsS3",
      "typeVersion": 1,
      "position": [1120, 400],
      "credentials": {
        "aws": {
          "id": "aws_credentials",
          "name": "AWS S3 Access"
        }
      }
    }
  ],
  "connections": {
    "Analysis Complete Webhook": {
      "main": [
        [
          {
            "node": "Parse Analysis Data",
            "type": "main",
            "index": 0
          }
        ]
      ]
    },
    "Parse Analysis Data": {
      "main": [
        [
          {
            "node": "Store to Database",
            "type": "main",
            "index": 0
          }
        ]
      ]
    },
    "Store to Database": {
      "main": [
        [
          {
            "node": "Check Alert Needed",
            "type": "main",
            "index": 0
          }
        ]
      ]
    },
    "Check Alert Needed": {
      "main": [
        [
          {
            "node": "Send Alert Email",
            "type": "main",
            "index": 0
          },
          {
            "node": "Archive to S3",
            "type": "main",
            "index": 0
          }
        ],
        [
          {
            "node": "Archive to S3",
            "type": "main",
            "index": 0
          }
        ]
      ]
    }
  },
  "active": true,
  "settings": {},
  "id": "amr-analysis-workflow"
}
```

## 2. Real-Time Monitoring Workflow

```json
{
  "name": "AMR Real-Time Monitor",
  "nodes": [
    {
      "parameters": {
        "rule": {
          "interval": [
            {
              "field": "minutes",
              "minutesInterval": 15
            }
          ]
        }
      },
      "name": "Every 15 Minutes",
      "type": "n8n-nodes-base.cron",
      "typeVersion": 1,
      "position": [240, 300]
    },
    {
      "parameters": {
        "operation": "executeQuery",
        "query": "SELECT \n  COUNT(*) as total_analyses,\n  COUNT(CASE WHEN high_priority_findings > 0 THEN 1 END) as high_priority_count,\n  AVG(quality_score) as avg_quality,\n  COUNT(CASE WHEN timestamp > NOW() - INTERVAL '1 hour' THEN 1 END) as recent_analyses\nFROM amr_analyses\nWHERE timestamp > NOW() - INTERVAL '24 hours'"
      },
      "name": "Query Recent Activity",
      "type": "n8n-nodes-base.postgres",
      "typeVersion": 1,
      "position": [460, 300],
      "credentials": {
        "postgres": {
          "id": "amr_database",
          "name": "AMR Results Database"
        }
      }
    },
    {
      "parameters": {
        "conditions": {
          "number": [
            {
              "value1": "={{$json.high_priority_count}}",
              "operation": "larger",
              "value2": 10
            }
          ]
        }
      },
      "name": "Check High Activity",
      "type": "n8n-nodes-base.if",
      "typeVersion": 1,
      "position": [680, 300]
    },
    {
      "parameters": {
        "channel": "#amr-surveillance",
        "text": "🔬 AMR Surveillance Update (Last 24h):\n• Total Analyses: {{$json.total_analyses}}\n• High-Priority Findings: {{$json.high_priority_count}}\n• Average Quality: {{Math.round($json.avg_quality * 100)}}%\n• Recent Activity: {{$json.recent_analyses}} analyses in last hour",
        "otherOptions": {}
      },
      "name": "Post to Slack",
      "type": "n8n-nodes-base.slack",
      "typeVersion": 1,
      "position": [900, 200],
      "credentials": {
        "slackApi": {
          "id": "slack_webhook",
          "name": "Slack AMR Channel"
        }
      }
    }
  ],
  "connections": {
    "Every 15 Minutes": {
      "main": [
        [
          {
            "node": "Query Recent Activity",
            "type": "main",
            "index": 0
          }
        ]
      ]
    },
    "Query Recent Activity": {
      "main": [
        [
          {
            "node": "Check High Activity",
            "type": "main",
            "index": 0
          }
        ]
      ]
    },
    "Check High Activity": {
      "main": [
        [
          {
            "node": "Post to Slack",
            "type": "main",
            "index": 0
          }
        ]
      ]
    }
  },
  "active": true,
  "settings": {},
  "id": "amr-monitoring-workflow"
}
```

## 3. Batch Processing Workflow

```json
{
  "name": "AMR Batch Processor",
  "nodes": [
    {
      "parameters": {
        "path": "amr-batch-process",
        "options": {}
      },
      "name": "Batch Request Webhook",
      "type": "n8n-nodes-base.webhook",
      "typeVersion": 1,
      "position": [240, 300],
      "webhookId": "amr-batch-process"
    },
    {
      "parameters": {
        "functionCode": "// Process batch analysis request\nconst request = items[0].json;\n\n// Validate batch request\nif (!request.genome_list || !Array.isArray(request.genome_list)) {\n  throw new Error('Invalid genome list provided');\n}\n\n// Split into manageable chunks\nconst chunkSize = 10;  // Process 10 genomes at a time\nconst chunks = [];\n\nfor (let i = 0; i < request.genome_list.length; i += chunkSize) {\n  chunks.push({\n    chunk_id: Math.floor(i / chunkSize) + 1,\n    total_chunks: Math.ceil(request.genome_list.length / chunkSize),\n    genomes: request.genome_list.slice(i, i + chunkSize),\n    user_email: request.user_email,\n    genes_file: request.genes_file,\n    batch_id: request.batch_id || 'batch_' + Date.now()\n  });\n}\n\nreturn chunks.map(chunk => ({json: chunk}));"
      },
      "name": "Split into Chunks",
      "type": "n8n-nodes-base.function",
      "typeVersion": 1,
      "position": [460, 300]
    },
    {
      "parameters": {
        "url": "http://localhost:8000/api/analyze",
        "sendHeaders": true,
        "headerParameters": {
          "parameters": [
            {
              "name": "Content-Type",
              "value": "application/json"
            }
          ]
        },
        "sendBody": true,
        "bodyParameters": {
          "parameters": [
            {
              "name": "genomes",
              "value": "={{$json.genomes}}"
            },
            {
              "name": "genes_file",
              "value": "={{$json.genes_file}}"
            },
            {
              "name": "user_email",
              "value": "={{$json.user_email}}"
            },
            {
              "name": "batch_mode",
              "value": "true"
            }
          ]
        },
        "options": {}
      },
      "name": "Process Chunk",
      "type": "n8n-nodes-base.httpRequest",
      "typeVersion": 3,
      "position": [680, 300]
    },
    {
      "parameters": {
        "operation": "insert",
        "table": "batch_processing_log",
        "columns": "batch_id,chunk_id,total_chunks,status,timestamp",
        "values": "={{$json.batch_id}},={{$json.chunk_id}},={{$json.total_chunks}},completed,={{new Date().toISOString()}}"
      },
      "name": "Log Progress",
      "type": "n8n-nodes-base.postgres",
      "typeVersion": 1,
      "position": [900, 300],
      "credentials": {
        "postgres": {
          "id": "amr_database",
          "name": "AMR Results Database"
        }
      }
    }
  ],
  "connections": {
    "Batch Request Webhook": {
      "main": [
        [
          {
            "node": "Split into Chunks",
            "type": "main",
            "index": 0
          }
        ]
      ]
    },
    "Split into Chunks": {
      "main": [
        [
          {
            "node": "Process Chunk",
            "type": "main",
            "index": 0
          }
        ]
      ]
    },
    "Process Chunk": {
      "main": [
        [
          {
            "node": "Log Progress",
            "type": "main",
            "index": 0
          }
        ]
      ]
    }
  },
  "active": true,
  "settings": {},
  "id": "amr-batch-workflow"
}
```

## Installation Instructions

1. **Install n8n**:
   ```bash
   npm install n8n -g
   # OR
   docker run -it --rm --name n8n -p 5678:5678 n8nio/n8n
   ```

2. **Import Workflows**:
   - Copy the JSON workflow definitions above
   - In n8n interface, go to "Workflows" > "Import from JSON"
   - Paste each workflow definition

3. **Configure Credentials**:
   - PostgreSQL/MySQL database connection
   - SMTP email settings
   - AWS S3 credentials (optional)
   - Slack webhook (optional)

4. **Test Integration**:
   ```python
   # Test the integration
   from src.utils.n8n_integration import integrate_with_pipeline
   
   test_results = {
       'genome_count': 5,
       'resistance_gene_count': 12,
       'high_priority_count': 3,
       'user_email': 'test@example.com'
   }
   
   integration_result = integrate_with_pipeline(test_results)
   print(integration_result)
   ```